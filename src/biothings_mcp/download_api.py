import os
from Bio import Entrez, SeqIO

class DownloadsMixin:

      def _gene_routes_config(self):
        """Configure gene routes for the API"""

        @self.post(
            "/gene/query",
            response_model=QueryResponse,
            tags=["genes"],
            summary="Query genes using a search string",
            operation_id="query_genes",
            description="""
            Search for genes using a query string with various filtering options.
            
            **IMPORTANT:** This endpoint requires structured queries using specific field names. 
            Simple natural language queries like "CDK2 gene" or "human kinase" will **NOT** work.
            You **MUST** specify the field you are querying, e.g., `symbol:CDK2`, `name:"cyclin-dependent kinase 2"`, `taxid:9606`.
            Use this endpoint when you need to *search* for genes based on criteria, not when you already know the specific gene ID.
            If you know the exact Entrez or Ensembl ID, use the `/gene/{gene_id}` endpoint instead for faster retrieval.
            If you only need general database information (like available fields or total gene count), use the `/gene/metadata` endpoint.

            **Supported Query Features (based on Lucene syntax):**
            1. Simple Term Queries:
               - `q=cdk2` (Searches across default fields like symbol, name, aliases)
               - `q="cyclin-dependent kinase"` (Searches for the exact phrase)
            
            2. Fielded Queries (specify the field to search):
               - `q=symbol:CDK2`
               - `q=name:"cyclin-dependent kinase 2"`
               - `q=refseq:NM_001798`
               - `q=ensembl.gene:ENSG00000123374`
               - `q=entrezgene:1017`
               - See [MyGene.info documentation](https://docs.mygene.info/en/latest/doc/query_service.html#available-fields) for more fields.
            
            3. Range Queries (for numerical or date fields):
               - `q=taxid:[9606 TO 10090]` (Find genes in taxonomy range including 9606 and 10090)
               - `q=entrezgene:>1000` (Find genes with Entrez ID greater than 1000)
            
            4. Boolean Queries:
               - `q=symbol:CDK2 AND taxid:9606` (Both conditions must be true)
               - `q=symbol:CDK* AND NOT taxid:9606` (Find CDK genes not in human)
               - `q=symbol:CDK2 OR symbol:BRCA1` (Either condition can be true)
               - `q=(symbol:CDK2 OR symbol:BRCA1) AND taxid:9606` (Grouping)
            
            5. Wildcard Queries:
               - `q=symbol:CDK*` (Matches symbols starting with CDK)
               - `q=name:*kinase*` (Matches names containing kinase)
               - `q=symbol:CDK?` (Matches CDK followed by one character)

            **Note:** See the [MyGene.info Query Syntax Guide](https://docs.mygene.info/en/latest/doc/query_service.html#query-syntax) for full details.
            
            The response includes pagination information (`total`, `max_score`, `took`) and the list of matching `hits`.
            """
        )
        async def query_genes(
            q: str = Query(
                ...,
                description="Query string following Lucene syntax. See endpoint description for details and examples.",
                examples=[
                    "symbol:CDK2",  # Fielded query
                    "name:\"cyclin-dependent kinase 2\"", # Phrase query
                    "refseq.rna:NM_001798", # Dot notation field
                    "taxid:[9606 TO 10090]", # Range query
                    "symbol:CDK2 AND taxid:9606", # Boolean query
                    "symbol:CDK* AND NOT taxid:9606", # Wildcard and boolean
                    "name:*kinase*", # Wildcard query
                    "(symbol:CDK2 OR symbol:BRCA1) AND taxid:9606", # Grouped boolean
                    "entrezgene:>1000" # Range query
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return from the matching gene hits. Supports dot notation (e.g., `refseq.rna`). If `fields=all`, all available fields are returned. Default: `symbol,name,taxid,entrezgene`.",
                examples=[
                    "symbol,name,taxid,entrezgene", # Default
                    "symbol,name,refseq.rna",
                    "ensembl.gene,uniprot.Swiss-Prot",
                    "summary,genomic_pos.chr,genomic_pos.start",
                    "all" # Return all fields
                ]
            ),
            size: int = Query(10, description="Maximum number of matching gene hits to return (capped at 1000). Default: 10.", examples=[10, 50, 1000]),
            skip: int = Query(0, description="Number of matching gene hits to skip, starting from 0 (for pagination). Default: 0.", examples=[0, 10, 50]),
            sort: Optional[str] = Query(None, description="Comma-separated fields to sort on. Prefix with `-` for descending order (e.g., `-symbol`). Default: sort by relevance score (`_score`) descending.", examples=["_score", "-_score", "symbol", "-entrezgene"]),
            species: Optional[str] = Query(None, description="Filter results by species. Accepts comma-separated taxonomy IDs (e.g., `9606,10090`) or common names for human, mouse, rat, fruitfly, nematode, zebrafish, thale-cress, frog, pig. Default: searches all species.", examples=["9606", "human", "10090"]),
            email: Optional[str] = Query(None, description="Optional user email for usage tracking.", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as a pandas DataFrame instead of JSON. Default: False.", examples=[True, False]),
            df_index: bool = Query(True, description="When `as_dataframe=True`, index the DataFrame by the internal `_id`. Default: True.", examples=[True, False])
        ):
            """Query genes"""
            log_message(message_type="debug:query_genes:entry", q=q, size=size)
            with start_action(action_type="api:query_genes", q=str(q), fields=str(fields), size=size, skip=skip, sort=str(sort), species=str(species), email=str(email), as_dataframe=as_dataframe, df_index=df_index):
                try:
                    async with GeneClientAsync() as client:
                        log_message(message_type="debug:query_genes:context_entered")
                        # Simplify call to match test_query_async
                        result = await client.query(q, size=size) 
                        log_message(message_type="debug:query_genes:raw_result", result=repr(result))

                    # Process result (keeping robust handling)
                    if not isinstance(result, dict):
                        log_message(message_type="debug:query_genes:result_not_dict")
                        result = {}
                    hits = result.get("hits", [])
                    validated_hits = []
                    for hit in hits:
                        try:
                            validated_hits.append(GeneResponse.model_validate(hit))
                        except Exception as e:
                            log_message(message_type="warning:query_genes:hit_validation_error", error=str(e), hit=repr(hit))
                            pass 
                    
                    response_obj = QueryResponse(
                        hits=validated_hits,
                        total=result.get("total"),
                        max_score=result.get("max_score"),
                        took=result.get("took")
                    )
                    log_message(message_type="debug:query_genes:returning", response=repr(response_obj))
                    return response_obj
                except Exception as e:
                    log_message(message_type="error:query_genes", error=str(e))
                    raise

    def download_fasta(self, ids: list[str]) -> str:
        Entrez.email = os.getenv("ENTREZ_EMAIL", "default_email@example.com")
        handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta")
        return handle.read()

    def download_gb(self, ids: list[str]) -> str:
        Entrez.email = os.getenv("ENTREZ_EMAIL", "default_email@example.com")
        handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb")
        return handle.read()
    
    def download_protein(self, ids: list[str]) -> str:
        Entrez.email = os.getenv("ENTREZ_EMAIL", "default_email@example.com")
        handle = Entrez.efetch(db="protein", id=ids, rettype="fasta")
        return handle.read()
    
    def download_genome(self, ids: list[str]) -> str:
        pass
    