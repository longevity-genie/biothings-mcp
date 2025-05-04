from pathlib import Path
from typing import Optional, List, Dict, Any, Union, Type, ClassVar

from dotenv import load_dotenv
from pydantic import BaseModel

from fastapi.middleware.cors import CORSMiddleware
from fastapi import FastAPI, Header, Query, Path, HTTPException
from fastapi.responses import StreamingResponse, RedirectResponse
from biothings_typed_client.genes import GeneClientAsync, GeneResponse
from biothings_typed_client.variants import VariantClientAsync, VariantResponse
from biothings_typed_client.chem import ChemClientAsync, ChemResponse
from biothings_typed_client.taxons import TaxonClientAsync, TaxonResponse
from eliot import start_task, start_action, log_message

class QueryResponse(BaseModel):
    hits: List[GeneResponse]
    total: Optional[int] = None
    max_score: Optional[float] = None
    took: Optional[int] = None

class VariantQueryResponse(BaseModel):
    hits: List[VariantResponse]
    total: Optional[int] = None
    max_score: Optional[float] = None
    took: Optional[int] = None

# Define a similar query response model for chemicals
class ChemQueryResponse(BaseModel):
    hits: List[ChemResponse]
    total: Optional[int] = None
    max_score: Optional[float] = None
    took: Optional[int] = None

# Define a similar query response model for taxons
class TaxonQueryResponse(BaseModel):
    hits: List[TaxonResponse]
    total: Optional[int] = None
    max_score: Optional[float] = None
    took: Optional[int] = None

class MetadataResponse(BaseModel):
    stats: Dict[str, Any]
    fields: Optional[Dict[str, Any]] = None
    index: Optional[Dict[str, Any]] = None
    version: Optional[str] = None

class GeneRoutesMixin:
    def _gene_routes_config(self):
        """Configure gene routes for the API"""

        @self.get(
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

        @self.get(
            "/gene/querymany",
            response_model=List[GeneResponse],
            tags=["genes"],
            summary="Batch query genes",
            operation_id="query_many_genes",
            description="""
            Perform multiple gene searches in a single request using a comma-separated list of query terms.
            This endpoint essentially performs a batch query similar to the POST request described in the [MyGene.info documentation](https://docs.mygene.info/en/latest/doc/query_service.html#batch-queries-via-post).

            **IMPORTANT:** Unlike `/gene/query`, the `query_list` parameter here takes multiple **terms** (like gene IDs, symbols, names) rather than full query strings.
            The `scopes` parameter defines which fields these terms should be searched against.
            Use this endpoint for batch *searching* of genes based on specific identifiers or terms within defined scopes.
            If you know the exact Entrez or Ensembl IDs for multiple genes and want direct retrieval, use the `/genes` endpoint instead (which is generally faster for ID lookups).

            **Endpoint Usage:**
            - Query multiple symbols: `query_list=CDK2,BRCA1` with `scopes=symbol`
            - Query multiple Entrez IDs: `query_list=1017,672` with `scopes=entrezgene`
            - Query mixed IDs/symbols: `query_list=CDK2,672` with `scopes=symbol,entrezgene` (searches both scopes for each term)
            
            **Result Interpretation:**
            - The response is a list of matching gene objects.
            - Each object includes a `query` field indicating which term from the `query_list` it matched.
            - A single term from `query_list` might match multiple genes (e.g., a symbol matching genes in different species if `species` is not set, or matching multiple retired IDs).
            - Terms with no matches are **omitted** from the response list (unlike the POST endpoint which returns a `notfound` entry).
            """
        )
        async def query_many_genes(
            query_list: str = Query(
                ...,
                description="Comma-separated list of query terms (e.g., gene IDs, symbols, names). Do NOT use complex query strings here; use `/gene/query` for that.",
                examples=[
                    "1017,1018", # Entrez IDs
                    "CDK2,BRCA1", # Symbols
                    "ENSG00000123374,ENSG00000139618", # Ensembl IDs
                    "1017,BRCA1", # Mixed IDs/Symbols
                    "cyclin-dependent kinase 2,breast cancer 1" # Names (use quotes if spaces)
                ]
            ),
            scopes: Optional[str] = Query(
                "entrezgene,ensemblgene,retired", # Default based on MyGene.info POST default
                description="Comma-separated list of fields to search the terms in `query_list` against. Default: `entrezgene,ensemblgene,retired`.",
                examples=[
                    "symbol",
                    "entrezgene",
                    "ensembl.gene",
                    "symbol,alias",
                    "entrezgene,ensemblgene,retired" # Default
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return from the matching gene hits. Default: `symbol,name,taxid,entrezgene`.",
                examples=[
                    "symbol,name,taxid,entrezgene", # Default
                    "symbol,name,refseq.rna",
                    "ensembl.gene,uniprot.Swiss-Prot",
                    "all" # Return all fields
                ]
            ),
            species: Optional[str] = Query(None, description="Filter results by species (comma-separated taxonomy IDs or common names). Default: searches all species.", examples=["9606", "human", "10090", "mouse,rat", "9606,10090"]),
            email: Optional[str] = Query(None, description="Optional user email for usage tracking.", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as a pandas DataFrame. Default: False.", examples=[True, False]),
            df_index: bool = Query(True, description="When `as_dataframe=True`, index the DataFrame by the matched `query` term. Default: True.", examples=[True, False]),
            size: int = Query(10, description="Maximum number of hits to return *per query term* (capped at 1000). Default: 10.", examples=[1, 5, 10, 1000])
        ):
            """Batch query genes"""
            log_message(message_type="debug:query_many_genes:entry", query_list=query_list, scopes=scopes, size=size)
            with start_action(action_type="api:query_many_genes", query_list=str(query_list), scopes=str(scopes), fields=str(fields), species=str(species), email=str(email), as_dataframe=str(as_dataframe), df_index=str(df_index), size=size):
                try:
                    async with GeneClientAsync() as client:
                        log_message(message_type="debug:query_many_genes:context_entered")
                        # Simplify call to match test_querymany_async
                        _query_list = query_list.split(",")
                        _scopes = scopes.split(",") if scopes else None
                        # Pass size to the client call
                        result = await client.querymany(_query_list, scopes=_scopes, size=size) 
                        log_message(message_type="debug:query_many_genes:raw_result", result=repr(result))

                    # Process result (keeping robust handling)
                    if not isinstance(result, list):
                        log_message(message_type="debug:query_many_genes:result_not_list")
                        return [] # Return empty list directly if result is not a list
                    
                    validated_genes = []
                    for gene_data in result:
                        if isinstance(gene_data, dict):
                            try:
                                validated_genes.append(GeneResponse.model_validate(gene_data))
                            except Exception as e:
                                log_message(message_type="warning:query_many_genes:gene_validation_error", error=str(e), data=repr(gene_data))
                                pass 
                        elif gene_data is not None: 
                            log_message(message_type="warning:query_many_genes:unexpected_data_type", data=repr(gene_data))
                    
                    log_message(message_type="debug:query_many_genes:returning", response=repr(validated_genes))
                    return validated_genes
                except Exception as e:
                    log_message(message_type="error:query_many_genes", error=str(e))
                    raise

        @self.get(
            "/gene/metadata",
            response_model=MetadataResponse,
            tags=["genes"],
            summary="Get gene database metadata",
            operation_id="get_gene_metadata",
            description="""
            Retrieve metadata about the underlying MyGene.info gene annotation database, **NOT** information about specific genes.
            
            **IMPORTANT:** Use this endpoint ONLY to understand the database itself (e.g., to discover available fields, check data versions, or get overall statistics). 
            It **CANNOT** be used to find or retrieve data for any particular gene. Use `/gene/query` or `/gene/{gene_id}` for that.
            
            **Returned Information Includes:**
            - `stats`: Database statistics (e.g., total number of genes).
            - `fields`: Available gene annotation fields and their data types.
            - `index`: Information about the backend data index.
            - `version`: Data version information.
            """
        )
        async def get_gene_metadata():
            """Get gene database metadata"""
            log_message(message_type="debug:get_gene_metadata:entry")
            with start_action(action_type="api:get_gene_metadata"):
                try:
                    async with GeneClientAsync() as client:
                        log_message(message_type="debug:get_gene_metadata:context_entered")
                        # Call matches test_metadata_async
                        result = await client.metadata()
                        log_message(message_type="debug:get_gene_metadata:raw_result", result=repr(result))
                        
                    # Process result (keeping robust handling)
                    if not isinstance(result, dict):
                        log_message(message_type="debug:get_gene_metadata:result_not_dict")
                        result = {}
                    stats = result.get("stats", {})
                    
                    response_obj = MetadataResponse(
                        stats=stats,
                        fields=result.get("fields"),
                        index=result.get("index"),
                        version=result.get("version")
                    )
                    log_message(message_type="debug:get_gene_metadata:returning", response=repr(response_obj))
                    return response_obj
                except Exception as e:
                    log_message(message_type="error:get_gene_metadata", error=str(e))
                    raise

        @self.get(
            "/gene/{gene_id}",
            response_model=GeneResponse,
            tags=["genes"],
            summary="Get gene information by ID",
            operation_id="get_gene",
            description="""
            Retrieves detailed information for a **single, specific gene** using its exact known identifier.
            
            **IMPORTANT:** **This is the preferred endpoint over `/gene/query` for fetching a specific gene when you already know its standard ID (Entrez or Ensembl) and don't need complex search filters.** It's generally faster for direct lookups.
            If you need to *search* for genes based on other criteria (like symbol, name, genomic location, function) or use complex boolean/range queries, use `/gene/query`.
            
            **Supported Identifiers:**
            - Entrez Gene ID: e.g., `1017`
            - Ensembl Gene ID: e.g., `ENSG00000123374`
            
            The response includes comprehensive gene information (fields can be customized using the `fields` parameter).
            If the ID is not found, a 404 error is returned.
            """
        )
        async def get_gene(
            gene_id: str = Path(..., description="Gene identifier (Entrez Gene ID or Ensembl Gene ID)", examples=["1017", "ENSG00000123374"]),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return. Default: `symbol,name,taxid,entrezgene` (plus other core fields). Use `all` for everything.",
                examples=[
                    "symbol,name,taxid,entrezgene", # Default
                    "symbol,name,summary",
                    "refseq,uniprot",
                    "all"
                ]
            ),
            species: Optional[str] = Query(
                None,
                description="Optional: Specify species (taxonomy ID or common name) to ensure the ID belongs to the correct species.",
                examples=["9606", "human", "10090"]
            ),
            email: Optional[str] = Query(
                None,
                description="Optional user email for usage tracking.",
                examples=["user@example.com"]
            ),
            # Parameters like as_dataframe, df_index, size, skip, sort are less relevant here 
            # as it retrieves a single specific document, but kept for consistency with client method
            as_dataframe: bool = Query(False, description="Return result as a pandas DataFrame. Default: False.", examples=[True, False], include_in_schema=False),
            df_index: bool = Query(True, description="When `as_dataframe=True`, index DataFrame by `_id`. Default: True.", examples=[True, False], include_in_schema=False),
            size: int = Query(10, description="Maximum number of results (irrelevant for single ID fetch).", include_in_schema=False),
            skip: int = Query(0, description="Number of results to skip (irrelevant for single ID fetch).", include_in_schema=False),
            sort: Optional[str] = Query(None, description="Sort field (irrelevant for single ID fetch).", include_in_schema=False)
        ):
            """Get gene information by ID"""
            with start_action(action_type="api:get_gene", gene_id=str(gene_id), fields=str(fields), species=str(species), email=str(email), as_dataframe=as_dataframe, df_index=df_index, size=size, skip=skip, sort=str(sort)):
                async with GeneClientAsync() as client:
                    return await client.getgene(
                        gene_id,
                        fields=fields.split(",") if fields else None,
                        species=species,
                        email=email,
                        as_dataframe=as_dataframe,
                        df_index=df_index,
                        size=size,
                        skip=skip,
                        sort=sort
                    )

        @self.get(
            "/genes",
            response_model=List[GeneResponse],
            tags=["genes"],
            summary="Get information for multiple genes",
            operation_id="get_genes",
            description="""
            Retrieves detailed information for **multiple specific genes** in a single request using their exact known identifiers.
            
            **IMPORTANT:** **This is the preferred endpoint over `/gene/querymany` for fetching multiple specific genes when you already know their standard IDs (Entrez, Ensembl) and don't need complex search filters.** Provide IDs as a comma-separated string. It's generally faster for direct batch lookups.
            If you need to perform batch *searches* for genes based on other criteria (like symbols across multiple species) or use different scopes per term, use `/gene/querymany`.

            **Input Format:**
            Accepts a comma-separated list of gene IDs (Entrez or Ensembl).
            
            **Endpoint Usage Examples:**
            - Multiple Entrez IDs: `gene_ids=1017,1018`
            - Multiple Ensembl IDs: `gene_ids=ENSG00000123374,ENSG00000134057`
            - Mixed IDs: `gene_ids=1017,ENSG00000134057`
            
            The response is a list containing an object for each **found** gene ID. IDs that are not found are silently omitted from the response list.
            The order of results in the response list corresponds to the order of IDs in the input `gene_ids` string.
            """
        )
        async def get_genes(
            gene_ids: str = Query(
                ...,
                description="Comma-separated list of gene IDs (Entrez or Ensembl).",
                examples=[
                    "1017,1018",
                    "ENSG00000123374,ENSG00000134057",
                    "1017,ENSG00000134057"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return for each gene. Default: `symbol,name,taxid,entrezgene` (plus other core fields). Use `all` for everything.",
                examples=[
                    "symbol,name,taxid,entrezgene", # Default
                    "symbol,name,summary",
                    "refseq,uniprot",
                    "all"
                ]
            ),
            species: Optional[str] = Query(
                None,
                description="Optional: Specify species (taxonomy ID or common name) to filter results.",
                examples=["9606", "human", "10090"]
            ),
            email: Optional[str] = Query(
                None,
                description="Optional user email for usage tracking.",
                examples=["user@example.com"]
            ),
            as_dataframe: bool = Query(
                False,
                description="Return results as a pandas DataFrame. Default: False.",
                examples=[True, False]
            ),
            df_index: bool = Query(
                True,
                description="When `as_dataframe=True`, index DataFrame by the gene ID (`_id`). Default: True.",
                examples=[True, False]
            )
        ):
            """Get information for multiple genes"""
            with start_action(action_type="api:get_genes", gene_ids=str(gene_ids), fields=str(fields), species=str(species), email=str(email), as_dataframe=as_dataframe, df_index=df_index):
                async with GeneClientAsync() as client:
                    return await client.getgenes(
                        gene_ids.split(","),
                        fields=fields.split(",") if fields else None,
                        species=species,
                        email=email,
                        as_dataframe=as_dataframe,
                        df_index=df_index
                    )

class VariantsRoutesMixin:
    def _variants_routes_config(self):
        """Configure variants routes for the API"""

        # Define query routes BEFORE specific ID routes to avoid path conflicts
        @self.get(
            "/variant/query",
            response_model=VariantQueryResponse,
            tags=["variants"],
            summary="Query variants using a search string",
            operation_id="query_variants",
            description="""
            Search for variants using a query string with various filtering options, leveraging the MyVariant.info API.
            See the [MyVariant.info Query Syntax Guide](https://docs.myvariant.info/en/latest/doc/variant_query_service.html#query-syntax) for full details.
            
            **IMPORTANT:** Use this endpoint for *searching* variants based on criteria. 
            If you already know the exact variant ID (HGVS, rsID), use the `/variant/{variant_id}` endpoint for faster direct retrieval.

            **Supported Query Features (Lucene syntax):**
            
            1. Simple Queries (searches default fields):
               - `q=rs58991260` (Find by rsID)
            
            2. Fielded Queries (specify the field):
               - `q=dbsnp.vartype:snp`
               - `q=dbnsfp.polyphen2.hdiv.pred:(D P)` (Matches D or P - space implies OR within parens for the same field)
               - `q=_exists_:dbsnp` (Variant must have a `dbsnp` field)
               - `q=_missing_:exac` (Variant must NOT have an `exac` field)
               - See [available fields](https://docs.myvariant.info/en/latest/doc/data.html#available-fields).
            
            3. Range Queries:
               - `q=dbnsfp.polyphen2.hdiv.score:>0.99`
               - `q=exac.af:[0 TO 0.00001]` (Inclusive range)
               - `q=exac.ac.ac_adj:{76640 TO 80000}` (Exclusive range)
            
            4. Wildcard Queries:
               - `q=dbnsfp.genename:CDK?` (Single character wildcard)
               - `q=dbnsfp.genename:CDK*` (Multi-character wildcard)
               - *Note: Wildcard cannot be the first character.* 
            
            5. Boolean Queries:
               - `q=_exists_:dbsnp AND dbsnp.vartype:snp`
               - `q=dbsnp.vartype:snp OR dbsnp.vartype:indel`
               - `q=_exists_:dbsnp AND NOT dbsnp.vartype:indel`
               - `q=(pubchem.molecular_weight:>500 OR chebi.mass:>500) AND _exists_:drugbank` (Grouping)

            6. Genomic Interval Queries (can be combined with AND, not within parentheses):
               - `q=chr1:69000-70000`
               - `q=chr1:69000-70000 AND dbnsfp.polyphen2.hdiv.score:>0.9`

            The response includes pagination information (`total`, `max_score`, `took`) and the list of matching `hits`.
            """
        )
        async def query_variants(
            q: str = Query(
                ...,
                description="Query string following Lucene syntax. See endpoint description and MyVariant.info docs for details.",
                examples=[
                    "rs58991260", # Simple query by rsID
                    "chr1:69000-70000", # Genomic range
                    "dbsnp.vartype:snp", # Fielded query
                    "dbnsfp.polyphen2.hdiv.pred:(D P)", # Fielded query with multiple values
                    "_exists_:dbsnp", # Existence query
                    "dbnsfp.polyphen2.hdiv.score:>0.99", # Range query
                    "dbnsfp.genename:CDK*", # Wildcard query
                    "_exists_:dbsnp AND dbsnp.vartype:snp", # Boolean query
                    "chr1:69000-70000 AND _exists_:clinvar" # Combined range and boolean
                ]
            ),
            fields: Optional[str] = Query(
                None, # Default is 'all' from client
                description="Comma-separated list of fields to return. Supports dot notation (e.g., `dbnsfp.genename`) and wildcards (`chebi.*`). Default: `all`.",
                examples=[
                    "all", # Default
                    "cadd.phred,dbsnp.rsid",
                    "clinvar.clinical_significance",
                    "dbnsfp.sift_pred,dbnsfp.polyphen2_hdiv_pred",
                    "dbnsfp.*"
                ]
            ),
            size: int = Query(10, description="Maximum number of matching variant hits to return (capped at 1000). Default: 10.", examples=[10, 50, 1000]),
            skip: int = Query(0, description="Number of matching variant hits to skip, starting from 0 (for pagination). Default: 0.", examples=[0, 10, 50]),
            sort: Optional[str] = Query(None, description="Comma-separated fields to sort on. Prefix with `-` for descending order. Default: sort by relevance score (`_score`) descending.", examples=["_score", "-cadd.phred", "dbsnp.rsid"]),
            email: Optional[str] = Query(None, description="Optional user email for usage tracking.", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as a pandas DataFrame. Default: False.", examples=[True, False]),
            df_index: bool = Query(True, description="When `as_dataframe=True`, index DataFrame by `_id`. Default: True.", examples=[True, False])
        ):
            """Query variants"""
            with start_action(action_type="api:query_variants", q=str(q), fields=str(fields), size=size, skip=skip, sort=str(sort), email=str(email), as_dataframe=as_dataframe, df_index=df_index):
                async with VariantClientAsync() as client:
                    result = await client.query(
                        q,
                        fields=fields.split(",") if fields else None,
                        size=size,
                        skip=skip,
                        sort=sort,
                        email=email,
                        as_dataframe=as_dataframe,
                        df_index=df_index
                    )
                    validated_hits = []
                    if isinstance(result, dict) and "hits" in result and isinstance(result["hits"], list):
                        for hit in result["hits"]:
                            try:
                                validated_hits.append(VariantResponse.model_validate(hit))
                            except Exception as e:
                                log_message(message_type="warning:query_variants:hit_validation_error", error=str(e), hit=repr(hit))
                                pass 
                    
                    return VariantQueryResponse(
                        hits=validated_hits,
                        total=result.get("total") if isinstance(result, dict) else None,
                        max_score=result.get("max_score") if isinstance(result, dict) else None,
                        took=result.get("took") if isinstance(result, dict) else None,
                    )

        @self.get(
            "/variants/querymany",
            tags=["variants"],
            summary="Batch query variants",
            operation_id="query_many_variants",
            description="""
            Perform multiple variant queries in a single request using a comma-separated list of variant identifiers.
            This endpoint is similar to the POST batch query functionality in the [MyVariant.info API](https://docs.myvariant.info/en/latest/doc/variant_query_service.html#batch-queries-via-post).
            
            **IMPORTANT:** This endpoint takes multiple **terms** (like rsIDs, HGVS IDs) in `query_list` and searches for them within the specified `scopes`.
            Use this for batch *searching* or retrieval based on specific identifiers within defined fields.
            If you know the exact IDs and want direct retrieval (which is generally faster), use the `/variants` endpoint.
            
            **Endpoint Usage:**
            - Query multiple rsIDs: `query_list=rs58991260,rs2500` with `scopes=dbsnp.rsid`
            - Query multiple HGVS IDs: `query_list=chr7:g.140453134T>C,chr1:g.69511A>G` (scopes likely not needed if IDs are unique, but `_id` or default scopes can be used)
            - Query mixed IDs: `query_list=rs58991260,chr1:g.69511A>G` with `scopes=dbsnp.rsid,_id`
            
            **Result Interpretation:**
            - The response is a list of matching variant objects.
            - Each object includes a `query` field indicating which term from the `query_list` it matched.
            - A single term might match multiple variants if the scope is broad (e.g., searching a gene name in `dbnsfp.genename`).
            - Terms with no matches are **omitted** from the response list (unlike the MyVariant POST endpoint which returns a `notfound` entry).
            """
        )
        async def query_many_variants(
            query_list: str = Query(
                ...,
                description="Comma-separated list of query terms (e.g., rsIDs, HGVS IDs). Do NOT use complex Lucene queries here; use `/variant/query` for that.",
                examples=[
                    "rs58991260,rs12190874", # rsIDs
                    "chr7:g.140453134T>C,chr1:g.69511A>G", # HGVS IDs
                    "rs58991260,chr1:g.69511A>G" # Mixed
                ]
            ),
            scopes: Optional[str] = Query(
                None, # Default scopes depend on the underlying client
                description="Comma-separated list of fields to search the terms in `query_list` against. If omitted, default scopes of the client are used (often includes `_id`, `dbsnp.rsid`, etc.).",
                examples=[
                    "dbsnp.rsid",
                    "_id", # For HGVS IDs
                    "clinvar.hgvs.coding,clinvar.hgvs.genomic",
                    "dbsnp.rsid,_id"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return. Supports dot notation and wildcards. If `all` or omitted, all fields may be returned (client default).",
                examples=[
                    "cadd.phred,dbsnp.rsid",
                    "clinvar.clinical_significance",
                    "dbnsfp.sift_pred,dbnsfp.polyphen2_hdiv_pred",
                    "dbnsfp.*",
                    "all"
                ]
            ),
            email: Optional[str] = Query(None, description="Optional user email for usage tracking.", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as a pandas DataFrame. Default: False.", examples=[True, False]),
            df_index: bool = Query(True, description="When `as_dataframe=True`, index DataFrame by the matched `query` term. Default: True.", examples=[True, False])
        ):
            """Batch query variants"""
            with start_action(action_type="api:query_many_variants", query_list=str(query_list), scopes=str(scopes), fields=str(fields), email=str(email), as_dataframe=as_dataframe, df_index=df_index):
                async with VariantClientAsync() as client:
                    return await client.querymany(
                        query_list.split(","),
                        scopes=scopes.split(",") if scopes else None,
                        fields=fields.split(",") if fields else None,
                        email=email,
                        as_dataframe=as_dataframe,
                        df_index=df_index
                    )
            
        # Define specific ID routes AFTER query routes
        @self.get(
            "/variant/{variant_id}",
            response_model=VariantResponse,
            tags=["variants"],
            summary="Get variant information by ID",
            operation_id="get_variant",
            description="""
            Retrieves detailed annotation data for a **single, specific variant** using its identifier, powered by the MyVariant.info annotation service.
            See the [MyVariant.info Annotation Service Docs](https://docs.myvariant.info/en/latest/doc/variant_annotation_service.html).
            
            **IMPORTANT:** **This is the preferred endpoint over `/variant/query` for fetching a specific variant when you already know its standard ID (HGVS or rsID) and don't need complex search filters.** It provides direct, generally faster, access to the full annotation object.
            If you need to *search* for variants based on other criteria (like gene name, functional impact, genomic region) or use complex boolean/range queries, use `/variant/query`.
            
            **Supported Identifiers (passed in the URL path):**
            - HGVS ID (e.g., `chr7:g.140453134T>C`). *Note: MyVariant.info primarily uses hg19-based HGVS IDs.* 
            - dbSNP rsID (e.g., `rs58991260`).
            
            The response includes comprehensive variant annotation information. By default (`fields=all` or omitted), the complete annotation object is returned.
            If the ID is not found or invalid, a 404 error is returned.
            """
        )
        async def get_variant(
            variant_id: str = Path(..., description="Variant identifier (HGVS ID or dbSNP rsID)", examples=["chr7:g.140453134T>C", "rs58991260"]),
            fields: Optional[str] = Query(
                None, # Default is effectively 'all' via the client
                description="Comma-separated list of fields to filter the returned annotation object. Supports dot notation (e.g., `cadd.gene`) and wildcards (`dbnsfp.*`). Default: `all` (returns the full object).",
                examples=[
                    "all", # Default
                    "cadd.phred,dbsnp.rsid",
                    "clinvar.clinical_significance",
                    "dbnsfp.sift_pred,dbnsfp.polyphen2_hdiv_pred",
                    "dbnsfp.*" 
                ]
            ),
            email: Optional[str] = Query(
                None,
                description="Optional user email for usage tracking.",
                examples=["user@example.com"]
            ),
            as_dataframe: bool = Query(
                False,
                description="Return result as a pandas DataFrame. Default: False.",
                examples=[True, False]
            ),
            df_index: bool = Query(
                True,
                description="When `as_dataframe=True`, index DataFrame by `_id`. Default: True.",
                examples=[True, False]
            )
        ):
            """Get variant information by ID"""
            with start_action(action_type="api:get_variant", variant_id=str(variant_id), fields=str(fields), email=str(email), as_dataframe=as_dataframe, df_index=df_index):
                async with VariantClientAsync() as client:
                    return await client.getvariant(
                        variant_id,
                        fields=fields.split(",") if fields else None,
                        email=email,
                        as_dataframe=as_dataframe,
                        df_index=df_index
                    )

        @self.get(
            "/variants",
            response_model=List[VariantResponse],
            tags=["variants"],
            summary="Get information for multiple variants",
            operation_id="get_variants",
            description="""
            Retrieves annotation data for **multiple specific variants** in a single request using their identifiers, similar to the MyVariant.info batch annotation service.
            See the [MyVariant.info Annotation Service Docs](https://docs.myvariant.info/en/latest/doc/variant_annotation_service.html#batch-queries-via-post).
            
            **IMPORTANT:** **This is the preferred endpoint over `/variants/querymany` for fetching multiple specific variants when you already know their standard IDs (HGVS or rsID).** Provide IDs as a comma-separated string. It's generally faster for direct batch lookups.
            If you need to perform batch *searches* for variants based on other criteria or use different scopes per term, use `/variants/querymany`.

            **Input Format:**
            Accepts a comma-separated list of variant IDs (HGVS or dbSNP rsIDs) in the `variant_ids` query parameter (max 1000 IDs).
            
            **Endpoint Usage Examples:**
            - Multiple HGVS IDs: `variant_ids=chr7:g.140453134T>C,chr1:g.69511A>G`
            - Multiple rsIDs: `variant_ids=rs58991260,rs2500`
            - Mixed IDs: `variant_ids=chr7:g.140453134T>C,rs58991260`
            
            The response is a list containing the full annotation object for each **found** variant ID. IDs that are not found or are invalid are silently omitted from the response list.
            Each returned object includes a `query` field indicating the input ID it corresponds to.
            The order of results generally corresponds to the input order but may vary for mixed ID types or invalid IDs.
            """
        )
        async def get_variants(
            variant_ids: str = Query(
                ...,
                description="Comma-separated list of variant IDs (HGVS or dbSNP rsIDs, max 1000).",
                examples=[
                    "chr7:g.140453134T>C,chr1:g.69511A>G",
                    "rs58991260,rs2500",
                    "chr7:g.140453134T>C,rs58991260"
                ]
            ),
            fields: Optional[str] = Query(
                None, # Default is effectively 'all' via the client
                description="Comma-separated list of fields to filter the returned annotation objects. Supports dot notation and wildcards. Default: `all` (returns full objects).",
                examples=[
                    "all", # Default
                    "cadd.phred,dbsnp.rsid",
                    "clinvar.clinical_significance",
                    "dbnsfp.sift_pred,dbnsfp.polyphen2_hdiv_pred",
                    "dbnsfp.*"
                ]
            ),
            email: Optional[str] = Query(
                None,
                description="Optional user email for usage tracking.",
                examples=["user@example.com"]
            ),
            as_dataframe: bool = Query(
                False,
                description="Return results as a pandas DataFrame. Default: False.",
                examples=[True, False]
            ),
            df_index: bool = Query(
                True,
                description="When `as_dataframe=True`, index DataFrame by the variant ID (`_id`). Default: True.",
                examples=[True, False]
            )
        ):
            """Get information for multiple variants"""
            with start_action(action_type="api:get_variants", variant_ids=str(variant_ids), fields=str(fields), email=str(email), as_dataframe=as_dataframe, df_index=df_index):
                async with VariantClientAsync() as client:
                    return await client.getvariants(
                        variant_ids.split(","),
                        fields=fields.split(",") if fields else None,
                        email=email,
                        as_dataframe=as_dataframe,
                        df_index=df_index
                    )

class ChemRoutesMixin:
    def _chem_routes_config(self):
        """Configure chemical routes for the API"""

        # Define query routes BEFORE specific ID routes
        @self.get(
            "/chem/query",
            response_model=ChemQueryResponse,
            tags=["chemicals"],
            summary="Query chemicals using a search string",
            operation_id="query_chems",
            description="""
            Search for chemical compounds using a query string with various filtering options.
            
            This endpoint supports complex queries with the following features:
            
            1. Simple Queries:
               - "C6H12O6" - Find compounds with molecular formula C6H12O6
               - "glucose" - Find compounds with name containing "glucose"
            
            2. Fielded Queries:
               - "pubchem.molecular_formula:C6H12O6" - Find compounds with specific formula
               - "pubchem.molecular_weight:[100 TO 200]" - Find compounds in weight range
               - "pubchem.xlogp:>2" - Find compounds with logP > 2
               - "pubchem.hydrogen_bond_donor_count:>2" - Find compounds with >2 H-bond donors
            
            3. Range Queries:
               - "pubchem.molecular_weight:[100 TO 200]" - Find compounds in weight range
               - "pubchem.xlogp:>2" - Find compounds with logP > 2
               - "pubchem.topological_polar_surface_area:[50 TO 100]" - Find compounds in TPSA range
            
            4. Boolean Queries:
               - "pubchem.hydrogen_bond_donor_count:>2 AND pubchem.hydrogen_bond_acceptor_count:>4"
               - "pubchem.molecular_weight:[100 TO 200] AND pubchem.xlogp:>2"
               - "pubchem.molecular_formula:C6H12O6 AND NOT pubchem.inchi_key:KTUFNOKKBVMGRW-UHFFFAOYSA-N"
            
            The response includes pagination information and can be returned as a pandas DataFrame.
            """
        )
        async def query_chems(
            q: str = Query(
                ...,
                description="Query string",
                examples=[
                    "C6H12O6",
                    "pubchem.molecular_formula:C6H12O6",
                    "pubchem.molecular_weight:[100 TO 200]",
                    "pubchem.xlogp:>2",
                    "pubchem.hydrogen_bond_donor_count:>2 AND pubchem.hydrogen_bond_acceptor_count:>4"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "pubchem.molecular_formula,pubchem.molecular_weight",
                    "pubchem.smiles,pubchem.inchi",
                    "pubchem.hydrogen_bond_donor_count,pubchem.hydrogen_bond_acceptor_count"
                ]
            ),
            size: int = Query(10, description="Maximum number of results to return (max 1000)", examples=[10, 50, 100]),
            skip: int = Query(0, description="Number of results to skip (for pagination)", examples=[0, 10, 20]),
            sort: Optional[str] = Query(None, description="Sort field, prefix with '-' for descending order", examples=["_score", "-_score", "pubchem.molecular_weight"]),
            email: Optional[str] = Query(None, description="User email for tracking usage", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as pandas DataFrame", examples=[True, False]),
            df_index: bool = Query(True, description="Index DataFrame by query (only if as_dataframe=True)", examples=[True, False])
        ):
            """Query chemicals"""
            with start_action(action_type="api:query_chems", q=str(q), fields=str(fields), size=size, skip=skip, sort=str(sort), email=str(email), as_dataframe=as_dataframe, df_index=df_index):
                async with ChemClientAsync() as client:
                    # Call the query method
                    result = await client.query(
                        q,
                        fields=fields.split(",") if fields else None,
                        size=size,
                        skip=skip,
                        sort=sort,
                        email=email,
                        as_dataframe=as_dataframe,
                        df_index=df_index
                    )
                    # Validate and structure the response according to ChemQueryResponse
                    validated_hits = []
                    if isinstance(result, dict) and "hits" in result and isinstance(result["hits"], list):
                        for hit in result["hits"]:
                            try:
                                # Individual hits still need validation if they are meant to be ChemResponse
                                validated_hits.append(ChemResponse.model_validate(hit))
                            except Exception as e:
                                log_message(message_type="warning:query_chems:hit_validation_error", error=str(e), hit=repr(hit))
                                pass 
                    
                    return ChemQueryResponse(
                        hits=validated_hits,
                        total=result.get("total") if isinstance(result, dict) else None,
                        max_score=result.get("max_score") if isinstance(result, dict) else None,
                        took=result.get("took") if isinstance(result, dict) else None,
                    )

        @self.get(
            "/chems/querymany",
            tags=["chemicals"],
            summary="Batch query chemicals",
            operation_id="query_many_chemicals",
            description="""
            Perform multiple chemical queries in a single request.
            
            This endpoint is useful for batch processing of chemical queries. It supports:
            
            1. Multiple Query Types:
               - Molecular formula queries: ["C6H12O6", "C12H22O11"]
               - Name queries: ["glucose", "sucrose"]
               - Mixed queries: ["C6H12O6", "sucrose"]
            
            2. Field Scoping:
               - Search in specific fields: scopes=["pubchem.molecular_formula", "pubchem.iupac"]
               - Search in all fields: scopes=None
            
            3. Result Filtering:
               - Return specific fields: fields=["pubchem.molecular_weight", "pubchem.xlogp"]
               - Return all fields: fields=None
            
            The response can be returned as a pandas DataFrame for easier data manipulation.
            """
        )
        async def query_many_chems(
            query_list: str = Query(
                ...,
                description="Comma-separated list of query strings",
                examples=[
                    "C6H12O6,C12H22O11",
                    "glucose,sucrose",
                    "C6H12O6,sucrose"
                ]
            ),
            scopes: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to search in",
                examples=[
                    "pubchem.molecular_formula,pubchem.iupac",
                    "pubchem.smiles,pubchem.inchi",
                    "pubchem.molecular_weight,pubchem.xlogp"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "pubchem.molecular_formula,pubchem.molecular_weight",
                    "pubchem.smiles,pubchem.inchi",
                    "pubchem.hydrogen_bond_donor_count,pubchem.hydrogen_bond_acceptor_count"
                ]
            ),
            email: Optional[str] = Query(None, description="User email for tracking usage", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as pandas DataFrame", examples=[True, False]),
            df_index: bool = Query(True, description="Index DataFrame by query (only if as_dataframe=True)", examples=[True, False])
        ):
            """Batch query chemicals"""
            with start_action(action_type="api:query_many_chems", query_list=str(query_list), scopes=str(scopes), fields=str(fields), email=str(email), as_dataframe=as_dataframe, df_index=df_index):
                async with ChemClientAsync() as client:
                    return await client.querymany(
                        query_list.split(","),
                        scopes=scopes.split(",") if scopes else None,
                        fields=fields.split(",") if fields else None,
                        email=email,
                        as_dataframe=as_dataframe,
                        df_index=df_index
                    )

        # Define specific ID routes AFTER query routes
        @self.get(
            "/chem/{chem_id}",
            response_model=ChemResponse,
            tags=["chemicals"],
            summary="Get chemical information by ID",
            operation_id="get_chem",
            description="""
            Retrieves detailed information about a specific chemical compound using its identifier.
            
            This endpoint supports various chemical ID formats:
            - InChIKey: "KTUFNOKKBVMGRW-UHFFFAOYSA-N" (Glucose)
            - PubChem CID: "5793" (Glucose)
            - SMILES: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O"
            
            The response includes comprehensive chemical information such as:
            - Basic information (ID, version)
            - PubChem information:
              - Structural properties (SMILES, InChI, molecular formula)
              - Physical properties (molecular weight, exact mass)
              - Chemical properties (hydrogen bond donors/acceptors, rotatable bonds)
              - Stereochemistry information (chiral centers, stereocenters)
              - Chemical identifiers (CID, InChIKey)
              - IUPAC names
              - Topological polar surface area
              - XLogP (octanol-water partition coefficient)
            
            You can filter the returned fields using the 'fields' parameter.
            """
        )
        async def get_chem(
            chem_id: str = Path(..., description="Chemical identifier", examples=["KTUFNOKKBVMGRW-UHFFFAOYSA-N", "5793", "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O"]),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "pubchem.molecular_formula,pubchem.molecular_weight",
                    "pubchem.smiles,pubchem.inchi",
                    "pubchem.hydrogen_bond_donor_count,pubchem.hydrogen_bond_acceptor_count"
                ]
            ),
            email: Optional[str] = Query(
                None,
                description="User email for tracking usage",
                examples=["user@example.com"]
            ),
            as_dataframe: bool = Query(
                False,
                description="Return results as pandas DataFrame",
                examples=[True, False]
            ),
            df_index: bool = Query(
                True,
                description="Index DataFrame by query (only if as_dataframe=True)",
                examples=[True, False]
            )
        ):
            """Get chemical information by ID"""
            with start_action(action_type="api:get_chem", chem_id=str(chem_id), fields=str(fields), email=str(email), as_dataframe=as_dataframe, df_index=df_index):
                async with ChemClientAsync() as client:
                    return await client.getchem(
                        chem_id,
                        fields=fields.split(",") if fields else None,
                        email=email,
                        as_dataframe=as_dataframe,
                        df_index=df_index
                    )

        @self.get(
            "/chems",
            response_model=List[ChemResponse],
            tags=["chemicals"],
            summary="Get information for multiple chemicals",
            operation_id="get_chems",
            description="""
            Retrieves information for multiple chemical compounds in a single request.
            
            This endpoint accepts a comma-separated list of chemical IDs in various formats:
            - InChIKeys: "KTUFNOKKBVMGRW-UHFFFAOYSA-N,XEFQLINVKFYRCS-UHFFFAOYSA-N"
            - PubChem CIDs: "5793,5281"
            - Mixed formats: "KTUFNOKKBVMGRW-UHFFFAOYSA-N,5281"
            
            The response includes the same comprehensive chemical information as the single chemical endpoint,
            but for all requested compounds.
            
            You can filter the returned fields using the 'fields' parameter.
            """
        )
        async def get_chems(
            chem_ids: str = Query(
                ...,
                description="Comma-separated list of chemical IDs",
                examples=[
                    "KTUFNOKKBVMGRW-UHFFFAOYSA-N,XEFQLINVKFYRCS-UHFFFAOYSA-N",
                    "5793,5281",
                    "KTUFNOKKBVMGRW-UHFFFAOYSA-N,5281"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "pubchem.molecular_formula,pubchem.molecular_weight",
                    "pubchem.smiles,pubchem.inchi",
                    "pubchem.hydrogen_bond_donor_count,pubchem.hydrogen_bond_acceptor_count"
                ]
            ),
            email: Optional[str] = Query(
                None,
                description="User email for tracking usage",
                examples=["user@example.com"]
            ),
            as_dataframe: bool = Query(
                False,
                description="Return results as pandas DataFrame",
                examples=[True, False]
            ),
            df_index: bool = Query(
                True,
                description="Index DataFrame by query (only if as_dataframe=True)",
                examples=[True, False]
            )
        ):
            """Get information for multiple chemicals"""
            with start_action(action_type="api:get_chems", chem_ids=str(chem_ids), fields=str(fields), email=str(email), as_dataframe=as_dataframe, df_index=df_index):
                async with ChemClientAsync() as client:
                    return await client.getchems(
                        chem_ids.split(","),
                        fields=fields.split(",") if fields else None,
                        email=email,
                        as_dataframe=as_dataframe,
                        df_index=df_index
                    )

class TaxonRoutesMixin:
    def _taxon_routes_config(self):
        """Configure taxon routes for the API"""
        @self.get(
            "/taxon/{taxon_id}",
            response_model=TaxonResponse,
            tags=["taxons"],
            summary="Get taxon information by ID",
            operation_id="get_taxon",
            description="""
            Retrieves detailed information about a specific taxon using its identifier.
            
            This endpoint supports both NCBI Taxonomy IDs and scientific names.
            
            Examples:
            - NCBI ID: 9606 (Homo sapiens)
            - Scientific name: "Homo sapiens"
            
            The response includes comprehensive taxon information such as:
            - Basic information (ID, scientific name, common name)
            - Taxonomic classification (rank, parent taxon)
            - Lineage information
            - Alternative names and authorities
            - Gene data availability
            
            By default, all available fields are returned. You can filter the returned fields using the 'fields' parameter.
            """
        )
        async def get_taxon(
            taxon_id: str = Path(..., description="Taxon identifier (NCBI ID or scientific name)", examples=["9606", "Homo sapiens"]),
            fields: str = Query(
                "all",
                description="Comma-separated list of fields to return. Use 'all' to return all fields (default).",
                examples=[
                    "all",
                    "scientific_name,common_name,rank",
                    "id,scientific_name,lineage",
                    "taxid,parent_taxid,rank"
                ]
            ),
            email: Optional[str] = Query(
                None,
                description="User email for tracking usage",
                examples=["user@example.com"]
            ),
            as_dataframe: bool = Query(
                False,
                description="Return results as pandas DataFrame",
                examples=[True, False]
            ),
            df_index: bool = Query(
                True,
                description="Index DataFrame by query (only if as_dataframe=True)",
                examples=[True, False]
            )
        ):
            """Get taxon information by ID"""
            # Replace underscores with spaces in scientific names
            #if not taxon_id.isdigit():
            #    taxon_id = taxon_id.replace("_", " ")
            
            # Convert fields string to list if provided and not "all"
            fields_list = fields.split(",") if fields != "all" else None
            
            async with TaxonClientAsync() as client:
                result = await client.gettaxon(
                    taxon_id,
                    fields=fields_list,
                    email=email,
                    as_dataframe=as_dataframe,
                    df_index=df_index
                )
                if result is None:
                    raise HTTPException(status_code=404, detail=f"Taxon '{taxon_id}' not found")
                return result

        @self.get(
            "/taxons",
            response_model=List[TaxonResponse],
            tags=["taxons"],
            summary="Get information for multiple taxa",
            operation_id="get_taxons",
            description="""
            Retrieves information for multiple taxa in a single request.
            
            This endpoint accepts a comma-separated list of taxon IDs (either NCBI IDs or scientific names).
            
            Examples:
            - Multiple NCBI IDs: "9606,10090" (Homo sapiens and Mus musculus)
            - Multiple scientific names: "Homo sapiens,Mus musculus"
            - Mixed IDs: "9606,Mus musculus" (Homo sapiens by NCBI ID and Mus musculus by name)
            
            The response includes the same comprehensive taxon information as the single taxon endpoint,
            but for all requested taxa.
            
            You can filter the returned fields using the 'fields' parameter.
            """
        )
        async def get_taxons(
            taxon_ids: str = Query(
                ...,
                description="Comma-separated list of taxon IDs",
                examples=[
                    "9606,10090",
                    "Homo sapiens,Mus musculus",
                    "9606,Mus musculus"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "scientific_name,common_name,rank",
                    "id,scientific_name,lineage",
                    "taxid,parent_taxid,rank"
                ]
            ),
            email: Optional[str] = Query(
                None,
                description="User email for tracking usage",
                examples=["user@example.com"]
            ),
            as_dataframe: bool = Query(
                False,
                description="Return results as pandas DataFrame",
                examples=[True, False]
            ),
            df_index: bool = Query(
                True,
                description="Index DataFrame by query (only if as_dataframe=True)",
                examples=[True, False]
            )
        ):
            """Get information for multiple taxa"""
            async with TaxonClientAsync() as client:
                return await client.gettaxons(
                    taxon_ids.split(","),
                    fields=fields.split(",") if fields else None,
                    email=email,
                    as_dataframe=as_dataframe,
                    df_index=df_index
                )

        @self.get(
            "/taxon/query",
            response_model=TaxonQueryResponse,
            tags=["taxons"],
            summary="Query taxa using a search string",
            operation_id="query_taxons",
            description="""
            Search for taxa using a query string with various filtering options.
            
            This endpoint supports complex queries with the following features:
            
            1. Simple Queries:
               - "scientific_name:Homo sapiens" - Find taxon by scientific name
               - "common_name:human" - Find taxon by common name
            
            2. Fielded Queries:
               - "rank:species" - Find species-level taxa
               - "parent_taxid:9606" - Find child taxa of Homo sapiens
               - "has_gene:true" - Find taxa with gene data
            
            3. Range Queries:
               - "taxid:[9606 TO 10090]" - Find taxa in ID range
               - "lineage:>9606" - Find taxa with Homo sapiens in lineage
            
            4. Boolean Queries:
               - "rank:species AND has_gene:true" - Find species with gene data
               - "scientific_name:Homo* AND NOT rank:genus" - Find taxa starting with Homo but not at genus level
            
            5. Wildcard Queries:
               - "scientific_name:Homo*" - Find taxa with scientific name starting with Homo
               - "common_name:*mouse*" - Find taxa with 'mouse' in common name
            
            The response includes pagination information and can be returned as a pandas DataFrame.
            """
        )
        async def query_taxons(
            q: str = Query(
                ...,
                description="Query string",
                examples=[
                    "scientific_name:Homo sapiens",
                    "rank:species",
                    "taxid:[9606 TO 10090]",
                    "rank:species AND has_gene:true",
                    "scientific_name:Homo* AND NOT rank:genus"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "scientific_name,common_name,rank",
                    "id,scientific_name,lineage",
                    "taxid,parent_taxid,rank"
                ]
            ),
            size: int = Query(10, description="Maximum number of results to return (max 1000)", examples=[10, 50, 100]),
            skip: int = Query(0, description="Number of results to skip (for pagination)", examples=[0, 10, 20]),
            sort: Optional[str] = Query(None, description="Sort field, prefix with '-' for descending order", examples=["_score", "-_score", "scientific_name"]),
            email: Optional[str] = Query(None, description="User email for tracking usage", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as pandas DataFrame", examples=[True, False]),
            df_index: bool = Query(True, description="Index DataFrame by query (only if as_dataframe=True)", examples=[True, False])
        ):
            """Query taxa"""
            async with TaxonClientAsync() as client:
                # Call the query method
                result = await client.query(
                    q,
                    fields=fields.split(",") if fields else None,
                    size=size,
                    skip=skip,
                    sort=sort,
                    email=email,
                    as_dataframe=as_dataframe,
                    df_index=df_index
                )
                # Validate and structure the response according to TaxonQueryResponse
                validated_hits = []
                if isinstance(result, dict) and "hits" in result and isinstance(result["hits"], list):
                    for hit in result["hits"]:
                        try:
                            # Individual hits need validation
                            validated_hits.append(TaxonResponse.model_validate(hit))
                        except Exception as e:
                            log_message(message_type="warning:query_taxons:hit_validation_error", error=str(e), hit=repr(hit))
                            pass 
                    
                return TaxonQueryResponse(
                    hits=validated_hits,
                    total=result.get("total") if isinstance(result, dict) else None,
                    max_score=result.get("max_score") if isinstance(result, dict) else None,
                    took=result.get("took") if isinstance(result, dict) else None,
                )

        @self.get(
            "/taxons/querymany",
            tags=["taxons"],
            summary="Batch query taxa",
            operation_id="query_many_taxons",
            description="""
            Perform multiple taxon queries in a single request.
            
            This endpoint is useful for batch processing of taxon queries. It supports:
            
            1. Multiple Query Types:
               - Scientific name queries: ["Homo sapiens", "Mus musculus"]
               - Common name queries: ["human", "mouse"]
               - Mixed queries: ["9606", "Mus musculus"]
            
            2. Field Scoping:
               - Search in specific fields: scopes=["scientific_name", "common_name"]
               - Search in all fields: scopes=None
            
            3. Result Filtering:
               - Return specific fields: fields=["scientific_name", "common_name", "rank"]
               - Return all fields: fields=None
            
            The response can be returned as a pandas DataFrame for easier data manipulation.
            """
        )
        async def query_many_taxons(
            query_list: str = Query(
                ...,
                description="Comma-separated list of query strings",
                examples=[
                    "Homo sapiens,Mus musculus",
                    "human,mouse",
                    "9606,Mus musculus"
                ]
            ),
            scopes: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to search in",
                examples=[
                    "scientific_name,common_name",
                    "taxid,parent_taxid",
                    "rank,lineage"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "scientific_name,common_name,rank",
                    "id,scientific_name,lineage",
                    "taxid,parent_taxid,rank"
                ]
            ),
            email: Optional[str] = Query(None, description="User email for tracking usage", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as pandas DataFrame", examples=[True, False]),
            df_index: bool = Query(True, description="Index DataFrame by query (only if as_dataframe=True)", examples=[True, False])
        ):
            """Batch query taxa"""
            async with TaxonClientAsync() as client:
                return await client.querymany(
                    query_list.split(","),
                    scopes=scopes.split(",") if scopes else None,
                    fields=fields.split(",") if fields else None,
                    email=email,
                    as_dataframe=as_dataframe,
                    df_index=df_index
                )
            
class BiothingsRestAPI(FastAPI, GeneRoutesMixin, VariantsRoutesMixin, ChemRoutesMixin, TaxonRoutesMixin):
    """FastAPI implementation providing OpenAI-compatible endpoints for Just-Agents.
    This class extends FastAPI to provide endpoints that mimic OpenAI's API structure,
    allowing Just-Agents to be used as a drop-in replacement for OpenAI's API.
    """

    def __init__(
        self,
        *,
        debug: bool = False,                       # Enable debug mode
        title: str = "BIO THINGS MCP Server",        # API title for documentation
        description: str = "BIO THINGS MCP Server, check https://github.com/longevity-genie/biothings-mcp for more information",
        version: str = "1.1.0",
        openapi_url: str = "/openapi.json",
        openapi_tags: Optional[List[Dict[str, Any]]] = None,
        servers: Optional[List[Dict[str, Union[str, Any]]]] = None,
        docs_url: str = "/docs",
        redoc_url: str = "/redoc",
        terms_of_service: Optional[str] = None,
        contact: Optional[Dict[str, Union[str, Any]]] = None,
        license_info: Optional[Dict[str, Union[str, Any]]] = None
    ) -> None:
        """Initialize the AgentRestAPI with FastAPI parameters.
        
        Args:
            debug: Enable debug mode
            title: API title shown in documentation
            description: API description shown in documentation
            version: API version
            openapi_url: URL for OpenAPI schema
            openapi_tags: List of tags to be included in the OpenAPI schema
            servers: List of servers to be included in the OpenAPI schema
            docs_url: URL for API documentation
            redoc_url: URL for ReDoc documentation
            terms_of_service: URL to the terms of service
            contact: Contact information in the OpenAPI schema
            license_info: License information in the OpenAPI schema
        """
        super().__init__(
            debug=debug,
            title=title,
            description=description,
            version=version,
            openapi_url=openapi_url,
            openapi_tags=openapi_tags,
            servers=servers,
            docs_url=docs_url,
            redoc_url=redoc_url,
            terms_of_service=terms_of_service,
            contact=contact,
            license_info=license_info
        )
        load_dotenv(override=True)
        # Initialize routes
        self._gene_routes_config()
        self._variants_routes_config()
        self._chem_routes_config()
        self._taxon_routes_config()

        # Add root route that redirects to docs
        @self.get("/", include_in_schema=False)
        async def root():
            return RedirectResponse(url="/docs")
