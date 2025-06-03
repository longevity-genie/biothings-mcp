import logging
from typing import Optional, List, Dict, Any, Union

# Setup logger for this module
logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG) # Optional: set level

from dotenv import load_dotenv
from pydantic import BaseModel, Field, ConfigDict

from fastapi import FastAPI, Path as FastApiPath, HTTPException
from fastapi.responses import RedirectResponse
from biothings_typed_client.genes import GeneClientAsync, GeneResponse
from biothings_typed_client.variants import VariantClientAsync, VariantResponse
from biothings_typed_client.chem import ChemClientAsync, ChemResponse
from biothings_typed_client.taxons import TaxonClientAsync, TaxonResponse as BaseClientTaxonResponse
from biothings_mcp.download_api import DownloadsMixin
from eliot import start_action

# Custom TaxonResponse model making version field optional
class TaxonResponse(BaseClientTaxonResponse):
    """Response model for taxon information with version field as optional"""
    model_config = ConfigDict(extra='allow')
    
    # Override parent's fields to set new serialization_alias and modify type (for version)
    # Also ensure validation_alias is present for input mapping
    id: str = Field(validation_alias='_id', serialization_alias='_id', description="Taxon identifier")
    version: Optional[int] = Field(default=1, validation_alias='_version', serialization_alias='_version', description="Version number")

# Request Body Models

class GeneQueryRequest(BaseModel):
    q: str = Field(..., description="Query string following Lucene syntax.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    size: int = Field(10, description="Maximum number of hits.")
    skip: int = Field(0, description="Number of hits to skip.")
    sort: Optional[str] = Field(None, description="Comma-separated fields to sort on.")
    species: Optional[str] = Field(None, description="Filter by species.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by _id.")

class GeneQueryManyRequest(BaseModel):
    query_list: str = Field(..., description="Comma-separated list of query terms.")
    scopes: Optional[str] = Field("entrezgene,ensemblgene,retired", description="Comma-separated list of fields to search against.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    species: Optional[str] = Field(None, description="Filter by species.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by matched query term.")
    size: int = Field(10, description="Maximum number of hits per query term.")

class GeneMetadataRequest(BaseModel):
    # This endpoint typically doesn't require a body,
    # but if it were to, it would be defined here.
    # For now, let's assume it might take an optional email for consistency.
    email: Optional[str] = Field(None, description="Optional user email for usage tracking.")
    pass


class GeneRequest(BaseModel):
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    species: Optional[str] = Field(None, description="Specify species.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.", json_schema_extra=False) # Kept from original, though less common for single GETs
    df_index: bool = Field(True, description="Index DataFrame by _id.", json_schema_extra=False)
    # size, skip, sort are usually not applicable for single ID fetch but were in original Query params
    size: int = Field(10, description="Max results (less relevant).", json_schema_extra=False)
    skip: int = Field(0, description="Skip results (less relevant).", json_schema_extra=False)
    sort: Optional[str] = Field(None, description="Sort field (less relevant).", json_schema_extra=False)


class GenesRequest(BaseModel):
    gene_ids: str = Field(..., description="Comma-separated list of gene IDs.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    species: Optional[str] = Field(None, description="Filter by species.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by gene ID.")

class VariantQueryRequest(BaseModel):
    q: str = Field(..., description="Query string following Lucene syntax.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    size: int = Field(10, description="Maximum number of hits.")
    skip: int = Field(0, description="Number of hits to skip.")
    sort: Optional[str] = Field(None, description="Comma-separated fields to sort on.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by _id.")

class VariantQueryManyRequest(BaseModel):
    query_list: str = Field(..., description="Comma-separated list of query terms.")
    scopes: Optional[str] = Field(None, description="Comma-separated list of fields to search against.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by matched query term.")

class VariantRequest(BaseModel):
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by _id.")

class VariantsRequest(BaseModel):
    variant_ids: str = Field(..., description="Comma-separated list of variant IDs.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by variant ID.")

class ChemQueryRequest(BaseModel):
    q: str = Field(..., description="Query string.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    size: int = Field(10, description="Maximum number of results.")
    skip: int = Field(0, description="Number of results to skip.")
    sort: Optional[str] = Field(None, description="Sort field.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by query.")

class ChemQueryManyRequest(BaseModel):
    query_list: str = Field(..., description="Comma-separated list of query strings.")
    scopes: Optional[str] = Field(None, description="Comma-separated list of fields to search in.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by query.")

class ChemRequest(BaseModel):
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by query.")

class ChemsRequest(BaseModel):
    chem_ids: str = Field(..., description="Comma-separated list of chemical IDs.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by query.")

class TaxonRequest(BaseModel):
    fields: str = Field("all", description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by query.")

class TaxonsRequest(BaseModel):
    taxon_ids: str = Field(..., description="Comma-separated list of taxon IDs.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by query.")

class TaxonQueryRequest(BaseModel):
    q: str = Field(..., description="Query string.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    size: int = Field(10, description="Maximum number of results.")
    skip: int = Field(0, description="Number of results to skip.")
    sort: Optional[str] = Field(None, description="Sort field.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by query.")

class TaxonQueryManyRequest(BaseModel):
    query_list: str = Field(..., description="Comma-separated list of query strings.")
    scopes: Optional[str] = Field(None, description="Comma-separated list of fields to search in.")
    fields: Optional[str] = Field(None, description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by query.")


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

        @self.post(
            "/gene/query",
            response_model=QueryResponse,
            tags=["genes"],
            summary="Search genes via Lucene query, returning gene details and query metadata.",
            operation_id="query_genes",
            response_description="Provides a list of `GeneResponse` objects (e.g., symbol, name, taxid) as 'hits', along with 'total' count, 'max_score', and 'took' time.",
            description="""
            Search for genes using a query string with various filtering options.
            
            **IMPORTANT:** This endpoint requires structured queries using specific field names. 
            Simple natural language queries like "CDK2 gene" or "human kinase" will **NOT** work.
            You **MUST** specify the field you are querying, e.g., `symbol:CDK2`, `name:"cyclin-dependent kinase 2"`, `taxid:9606`.
            Use this endpoint when you need to *search* for genes based on criteria, not when you already know the specific gene ID.
            If you know the exact Entrez or Ensembl ID, use the `/gene/{gene_id}` endpoint instead for faster retrieval.
            If you only need general database information (like available fields or total gene count), use the `/gene/metadata` endpoint.
            It does not give exact sequences of the gene but gives ids of genomic, protein and rna sequences which you can download with other tools (like download_entrez_data). If you use those ids in downloads you must alway check whether user wants protein, dna or rna sequence (clarify if not clear)
            It does not give exact variants but you have tool for variants resolution based on this gene ids

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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def query_genes(request: GeneQueryRequest):
            """Query genes"""
            with start_action(action_type="api:query_genes", q=str(request.q), size=request.size) as action:
                try:
                    result = {}
                    async with GeneClientAsync() as client:
                        action.log(message_type="debug:query_genes:context_entered")
                        result = await client.query(
                            request.q,
                            fields=request.fields.split(",") if request.fields else None,
                            size=request.size,
                            skip=request.skip,
                            sort=request.sort,
                            species=request.species,
                            email=request.email,
                            as_dataframe=request.as_dataframe,
                            df_index=request.df_index
                        ) 
                        action.log(message_type="debug:query_genes:raw_result", result=repr(result))

                    if not isinstance(result, dict):
                        action.log(message_type="debug:query_genes:result_not_dict")
                        result = {}
                    hits = result.get("hits", [])
                    validated_hits = []
                    for hit in hits:
                        try:
                            validated_hits.append(GeneResponse.model_validate(hit))
                        except Exception as e:
                            action.log(message_type="warning:query_genes:hit_validation_error", error=str(e), hit=repr(hit))
                            pass 
                    
                    response_obj = QueryResponse(
                        hits=validated_hits,
                        total=result.get("total"),
                        max_score=result.get("max_score"),
                        took=result.get("took")
                    )
                    action.log(message_type="debug:query_genes:returning", response=repr(response_obj))
                    return response_obj
                except Exception as e:
                    action.log(message_type="error:query_genes", error=str(e))
                    raise

        @self.post(
            "/gene/querymany",
            response_model=List[GeneResponse],
            tags=["genes"],
            summary="Batch query genes by multiple terms, returning a list of gene details.",
            operation_id="query_many_genes",
            response_description="Returns a list of `GeneResponse` objects, each including details like symbol, name, taxid, and the original query term.",
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
            - It does not give exact sequences of the gene but gives ids of genomic, protein and rna sequences which you can download with other tools (like download_entrez_data). If you use those ids in downloads you must alway check whether user wants protein, dna or rna sequence (clarify if not clear)
            - It does not give exact variants but you have tool for variants resolution based on this gene ids
            - Each object includes a `query` field indicating which term from the `query_list` it matched.
            - A single term from `query_list` might match multiple genes (e.g., a symbol matching genes in different species if `species` is not set, or matching multiple retired IDs).
            - Terms with no matches are **omitted** from the response list (unlike the POST endpoint which returns a `notfound` entry).
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def query_many_genes(request: GeneQueryManyRequest):
            """Batch query genes"""
            with start_action(action_type="api:query_many_genes", query_list=str(request.query_list), scopes=str(request.scopes), fields=str(request.fields), species=str(request.species), email=str(request.email), as_dataframe=str(request.as_dataframe), df_index=str(request.df_index), size=request.size) as action:
                try:
                    result = []
                    async with GeneClientAsync() as client:
                        action.log(message_type="debug:query_many_genes:context_entered")
                        result = await client.querymany(
                            request.query_list.split(","), 
                            scopes=request.scopes.split(",") if request.scopes else None,
                            fields=request.fields.split(",") if request.fields else None,
                            species=request.species,
                            email=request.email,
                            as_dataframe=request.as_dataframe,
                            df_index=request.df_index,
                            size=request.size
                        ) 
                        action.log(message_type="debug:query_many_genes:raw_result", result=repr(result))

                    if not isinstance(result, list):
                        action.log(message_type="debug:query_many_genes:result_not_list")
                        return [] 
                    
                    validated_genes = []
                    for gene_data in result:
                        if isinstance(gene_data, dict):
                            try:
                                validated_genes.append(GeneResponse.model_validate(gene_data))
                            except Exception as e:
                                action.log(message_type="warning:query_many_genes:gene_validation_error", error=str(e), data=repr(gene_data))
                                pass 
                        elif gene_data is not None: 
                            action.log(message_type="warning:query_many_genes:unexpected_data_type", data=repr(gene_data))
                    
                    action.log(message_type="debug:query_many_genes:returning", response=repr(validated_genes))
                    return validated_genes
                except Exception as e:
                    action.log(message_type="error:query_many_genes", error=str(e))
                    raise

        @self.post(
            "/gene/metadata",
            response_model=MetadataResponse,
            tags=["genes"],
            summary="Retrieve MyGene.info database metadata including stats and fields.",
            operation_id="get_gene_metadata",
            response_description="Returns `MetadataResponse` containing database `stats` (e.g., total genes), available `fields` with data types, `index` information, and data `version`.",
            description="""
            Retrieve metadata about the underlying MyGene.info gene annotation database, **NOT** information about specific genes.
            
            **IMPORTANT:** Use this endpoint ONLY to understand the database itself (e.g., to discover available fields, check data versions, or get overall statistics). 
            It **CANNOT** be used to find or retrieve data for any particular gene. Use `/gene/query` or `/gene/{gene_id}` for that.
            
            **Returned Information Includes:**
            - `stats`: Database statistics (e.g., total number of genes).
            - `fields`: Available gene annotation fields and their data types.
            - `index`: Information about the backend data index.
            - `version`: Data version information.
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def get_gene_metadata(request: Optional[GeneMetadataRequest] = None): # Made request body optional
            """Get gene database metadata"""
            email_for_tracking = request.email if request and request.email else None
            # Log the email if provided, but don't pass it to client.metadata()
            with start_action(action_type="api:get_gene_metadata", email_provided=bool(email_for_tracking)) as action:
                try:
                    result = None
                    async with GeneClientAsync() as client:
                        action.log(message_type="debug:get_gene_metadata:context_entered")
                        # Call client.metadata() without the email argument
                        result = await client.metadata()
                        action.log(message_type="debug:get_gene_metadata:raw_result", result=repr(result))
                    
                    if not isinstance(result, dict):
                        action.log(message_type="debug:get_gene_metadata:result_not_dict")
                        result = {}
                    stats = result.get("stats", {})
                    
                    response_obj = MetadataResponse(
                        stats=stats,
                        fields=result.get("fields"),
                        index=result.get("index"),
                        version=result.get("version")
                    )
                    action.log(message_type="debug:get_gene_metadata:returning", response=repr(response_obj))
                    return response_obj
                except Exception as e:
                    action.log(message_type="error:get_gene_metadata", error=str(e))
                    raise

        @self.post(
            "/gene/{gene_id}",
            response_model=GeneResponse,
            tags=["genes"],
            summary="Fetch a specific gene by Entrez or Ensembl ID.",
            operation_id="get_gene",
            response_description="Returns a `GeneResponse` object containing detailed information for the specified gene, such as symbol, name, taxid, and entrezgene.",
            description="""
            Retrieves detailed information for a **single, specific gene** using its exact known identifier.
            
            **IMPORTANT:** **This is the preferred endpoint over `/gene/query` for fetching a specific gene when you already know its standard ID (Entrez or Ensembl) and don't need complex search filters.** It's generally faster for direct lookups.
            If you need to *search* for genes based on other criteria (like symbol, name, genomic location, function) or use complex boolean/range queries, use `/gene/query`.
            
            **Supported Identifiers:**
            - Entrez Gene ID: e.g., `1017`
            - Ensembl Gene ID: e.g., `ENSG00000123374`
            
            The response includes comprehensive gene information (fields can be customized using the `fields` parameter).
            If the ID is not found, a 404 error is returned.
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def get_gene(
            gene_id: str = FastApiPath(..., description="Gene identifier (Entrez Gene ID or Ensembl Gene ID)", examples=["1017", "ENSG00000123374"]),
            request: Optional[GeneRequest] = None # Made request body optional
        ):
            """Get gene information by ID"""
            # Use request body if provided, otherwise use defaults (empty or None)
            _fields = request.fields if request and request.fields else None
            _species = request.species if request and request.species else None
            _email = request.email if request and request.email else None
            _as_dataframe = request.as_dataframe if request else False
            _df_index = request.df_index if request else True
            _size = request.size if request else 10
            _skip = request.skip if request else 0
            _sort = request.sort if request and request.sort else None            

            with start_action(action_type="api:get_gene", gene_id=str(gene_id), fields=str(_fields), species=str(_species), email=str(_email), as_dataframe=_as_dataframe, df_index=_df_index, size=_size, skip=_skip, sort=str(_sort)) as action:
                async with GeneClientAsync() as client:
                    return await client.getgene(
                        gene_id,
                        fields=_fields.split(",") if _fields else None,
                        species=_species,
                        email=_email,
                        as_dataframe=_as_dataframe,
                        df_index=_df_index,
                        size=_size,
                        skip=_skip,
                        sort=_sort
                    )

        @self.post(
            "/genes",
            response_model=List[GeneResponse],
            tags=["genes"],
            summary="Fetch multiple genes by a comma-separated list of Entrez or Ensembl IDs.",
            operation_id="get_genes",
            response_description="Returns a list of `GeneResponse` objects, each containing detailed information (e.g., symbol, name, taxid) for the corresponding requested gene ID.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def get_genes(request: GenesRequest):
            """Get information for multiple genes"""
            with start_action(action_type="api:get_genes", gene_ids=str(request.gene_ids), fields=str(request.fields), species=str(request.species), email=str(request.email), as_dataframe=request.as_dataframe, df_index=request.df_index) as action:
                async with GeneClientAsync() as client:
                    return await client.getgenes(
                        request.gene_ids.split(","),
                        fields=request.fields.split(",") if request.fields else None,
                        species=request.species,
                        email=request.email,
                        as_dataframe=request.as_dataframe,
                        df_index=request.df_index
                    )

class VariantsRoutesMixin:
    def _variants_routes_config(self):
        """Configure variants routes for the API"""

        # Define query routes BEFORE specific ID routes to avoid path conflicts
        @self.post(
            "/variant/query",
            response_model=VariantQueryResponse,
            tags=["variants"],
            summary="Search variants via Lucene query (e.g., rsID, gene name), returning variant details and query metadata.",
            operation_id="query_variants",
            response_description="Provides a list of `VariantResponse` objects (e.g., id, chrom, vcf details) as 'hits', along with 'total' count, 'max_score', and 'took' time.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def query_variants(request: VariantQueryRequest):
            """Query variants"""
            with start_action(action_type="api:query_variants", q=str(request.q), fields=str(request.fields), size=request.size, skip=request.skip, sort=str(request.sort), email=str(request.email), as_dataframe=request.as_dataframe, df_index=request.df_index) as action:
                try:
                    result = {}
                    async with VariantClientAsync() as client:
                        result = await client.query(
                            request.q,
                            fields=request.fields.split(",") if request.fields else None,
                            size=request.size,
                            skip=request.skip,
                            sort=request.sort,
                            email=request.email,
                            as_dataframe=request.as_dataframe,
                            df_index=request.df_index
                        )
                        action.log(message_type="debug:query_variants:raw_result", result=repr(result))

                    validated_hits = []
                    if isinstance(result, dict) and "hits" in result and isinstance(result["hits"], list):
                        for hit in result["hits"]:
                            try:
                                validated_hits.append(VariantResponse.model_validate(hit))
                            except Exception as e:
                                action.log(message_type="warning:query_variants:hit_validation_error", error=str(e), hit=repr(hit))
                                pass 
                    
                    return VariantQueryResponse(
                        hits=validated_hits,
                        total=result.get("total") if isinstance(result, dict) else None,
                        max_score=result.get("max_score") if isinstance(result, dict) else None,
                        took=result.get("took") if isinstance(result, dict) else None,
                    )
                except Exception as e:
                    action.log(message_type="error:query_variants", error=str(e))
                    raise

        @self.post(
            "/variants/querymany",
            response_model=List[VariantResponse], # Assuming this based on client and pattern
            tags=["variants"],
            summary="Batch query variants by multiple identifiers (e.g., rsIDs, HGVS IDs).",
            operation_id="query_many_variants",
            response_description="Returns a list of `VariantResponse` objects, each including details like id, chrom, vcf information, and the original query term.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def query_many_variants(request: VariantQueryManyRequest):
            """Batch query variants"""
            with start_action(action_type="api:query_many_variants", query_list=str(request.query_list), scopes=str(request.scopes), fields=str(request.fields), email=str(request.email), as_dataframe=request.as_dataframe, df_index=request.df_index) as action:
                try:
                    result = []
                    async with VariantClientAsync() as client:
                        result = await client.querymany(
                            request.query_list.split(","),
                            scopes=request.scopes.split(",") if request.scopes else None,
                            fields=request.fields.split(",") if request.fields else None,
                            email=request.email,
                            as_dataframe=request.as_dataframe,
                            df_index=request.df_index
                        )
                        action.log(message_type="debug:query_many_variants:raw_result", result=repr(result))
                    return result
                except Exception as e:
                    action.log(message_type="error:query_many_variants", error=str(e))
                    raise
            
        # Define specific ID routes AFTER query routes
        @self.post(
            "/variant/{variant_id}",
            response_model=VariantResponse,
            tags=["variants"],
            summary="Fetch a specific variant by HGVS or rsID.",
            operation_id="get_variant",
            response_description="Returns a `VariantResponse` object containing detailed annotation data for the specified variant, such as id, chromosome, VCF info.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def get_variant(
            variant_id: str = FastApiPath(..., description="Variant identifier (HGVS ID or dbSNP rsID)", examples=["chr7:g.140453134T>C", "rs58991260"]),
            request: Optional[VariantRequest] = None # Made request body optional
        ):
            """Get variant information by ID"""
            _fields = request.fields if request and request.fields else None
            _email = request.email if request and request.email else None
            _as_dataframe = request.as_dataframe if request else False
            _df_index = request.df_index if request else True

            with start_action(action_type="api:get_variant", variant_id=str(variant_id), fields=str(_fields), email=str(_email), as_dataframe=_as_dataframe, df_index=_df_index) as action:
                async with VariantClientAsync() as client:
                    return await client.getvariant(
                        variant_id,
                        fields=_fields.split(",") if _fields else None,
                        email=_email,
                        as_dataframe=_as_dataframe,
                        df_index=_df_index
                    )

        @self.post(
            "/variants",
            response_model=List[VariantResponse],
            tags=["variants"],
            summary="Fetch multiple variants by a comma-separated list of HGVS or rsIDs.",
            operation_id="get_variants",
            response_description="Returns a list of `VariantResponse` objects, each containing detailed annotation data (e.g., id, chrom, vcf info) for the corresponding requested variant ID.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def get_variants(request: VariantsRequest):
            """Get information for multiple variants"""
            with start_action(action_type="api:get_variants", variant_ids=str(request.variant_ids), fields=str(request.fields), email=str(request.email), as_dataframe=request.as_dataframe, df_index=request.df_index) as action:
                async with VariantClientAsync() as client:
                    return await client.getvariants(
                        request.variant_ids.split(","),
                        fields=request.fields.split(",") if request.fields else None,
                        email=request.email,
                        as_dataframe=request.as_dataframe,
                        df_index=request.df_index
                    )

class ChemRoutesMixin:
    def _chem_routes_config(self):
        """Configure chemical routes for the API"""

        # Define query routes BEFORE specific ID routes
        @self.post(
            "/chem/query",
            response_model=ChemQueryResponse,
            tags=["chemicals"],
            summary="Search chemical compounds via Lucene query (e.g., name, formula), returning compound details and query metadata.",
            operation_id="query_chems",
            response_description="Provides a list of `ChemResponse` objects (e.g., id, PubChem formula, weight) as 'hits', along with 'total' count, 'max_score', and 'took' time.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def query_chems(request: ChemQueryRequest):
            """Query chemicals"""
            with start_action(action_type="api:query_chems", q=str(request.q), fields=str(request.fields), size=request.size, skip=request.skip, sort=str(request.sort), email=str(request.email), as_dataframe=request.as_dataframe, df_index=request.df_index) as action:
                async with ChemClientAsync() as client:
                    result = await client.query(
                        request.q,
                        fields=request.fields.split(",") if request.fields else None,
                        size=request.size,
                        skip=request.skip,
                        sort=request.sort,
                        email=request.email,
                        as_dataframe=request.as_dataframe,
                        df_index=request.df_index
                    )
                    validated_hits = []
                    if isinstance(result, dict) and "hits" in result and isinstance(result["hits"], list):
                        for hit in result["hits"]:
                            try:
                                validated_hits.append(ChemResponse.model_validate(hit))
                            except Exception as e:
                                action.log(message_type="warning:query_chems:hit_validation_error", error=str(e), hit=repr(hit))
                                pass 
                    
                    return ChemQueryResponse(
                        hits=validated_hits,
                        total=result.get("total") if isinstance(result, dict) else None,
                        max_score=result.get("max_score") if isinstance(result, dict) else None,
                        took=result.get("took") if isinstance(result, dict) else None,
                    )

        @self.post(
            "/chems/querymany",
            response_model=List[ChemResponse], # Assuming this based on client and pattern
            tags=["chemicals"],
            summary="Batch query chemical compounds by multiple terms (e.g., names, InChIKeys).",
            operation_id="query_many_chemicals",
            response_description="Returns a list of `ChemResponse` objects, each including details like id, PubChem information, and the original query term.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def query_many_chems(request: ChemQueryManyRequest):
            """Batch query chemicals"""
            with start_action(action_type="api:query_many_chems", query_list=str(request.query_list), scopes=str(request.scopes), fields=str(request.fields), email=str(request.email), as_dataframe=request.as_dataframe, df_index=request.df_index) as action:
                async with ChemClientAsync() as client:
                    return await client.querymany(
                        request.query_list.split(","),
                        scopes=request.scopes.split(",") if request.scopes else None,
                        fields=request.fields.split(",") if request.fields else None,
                        email=request.email,
                        as_dataframe=request.as_dataframe,
                        df_index=request.df_index
                    )

        # Define specific ID routes AFTER query routes
        @self.post(
            "/chem/{chem_id}",
            response_model=ChemResponse,
            tags=["chemicals"],
            summary="Fetch a specific chemical compound by ID (e.g., InChIKey, PubChem CID).",
            operation_id="get_chem",
            response_description="Returns a `ChemResponse` object containing detailed information for the specified chemical, including PubChem data like formula, weight, and XLogP.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def get_chem(
            chem_id: str = FastApiPath(..., description="Chemical identifier", examples=["KTUFNOKKBVMGRW-UHFFFAOYSA-N", "5793", "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O"]),
            request: Optional[ChemRequest] = None # Made request body optional
        ):
            """Get chemical information by ID"""
            _fields = request.fields if request and request.fields else None
            _email = request.email if request and request.email else None
            _as_dataframe = request.as_dataframe if request else False
            _df_index = request.df_index if request else True

            with start_action(action_type="api:get_chem", chem_id=str(chem_id), fields=str(_fields), email=str(_email), as_dataframe=_as_dataframe, df_index=_df_index) as action:
                async with ChemClientAsync() as client:
                    return await client.getchem(
                        chem_id,
                        fields=_fields.split(",") if _fields else None,
                        email=_email,
                        as_dataframe=_as_dataframe,
                        df_index=_df_index
                    )

        @self.post(
            "/chems",
            response_model=List[ChemResponse],
            tags=["chemicals"],
            summary="Fetch multiple chemical compounds by a comma-separated list of IDs.",
            operation_id="get_chems",
            response_description="Returns a list of `ChemResponse` objects, each containing detailed information (e.g., id, PubChem data) for the corresponding requested chemical ID.",
            description="""
            Retrieves information for multiple chemical compounds in a single request.
            
            This endpoint accepts a comma-separated list of chemical IDs in various formats:
            - InChIKeys: "KTUFNOKKBVMGRW-UHFFFAOYSA-N,XEFQLINVKFYRCS-UHFFFAOYSA-N"
            - PubChem CIDs: "5793,5281"
            - Mixed formats: "KTUFNOKKBVMGRW-UHFFFAOYSA-N,5281"
            
            The response includes the same comprehensive chemical information as the single chemical endpoint,
            but for all requested compounds.
            
            You can filter the returned fields using the 'fields' parameter.
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def get_chems(request: ChemsRequest):
            """Get information for multiple chemicals"""
            with start_action(action_type="api:get_chems", chem_ids=str(request.chem_ids), fields=str(request.fields), email=str(request.email), as_dataframe=request.as_dataframe, df_index=request.df_index) as action:
                async with ChemClientAsync() as client:
                    return await client.getchems(
                        request.chem_ids.split(","),
                        fields=request.fields.split(",") if request.fields else None,
                        email=request.email,
                        as_dataframe=request.as_dataframe,
                        df_index=request.df_index
                    )

class TaxonRoutesMixin:
    def _taxon_routes_config(self):
        """Configure taxon routes for the API"""
        @self.post(
            "/taxon/{taxon_id}",
            response_model=TaxonResponse,
            tags=["taxons"],
            summary="Fetch a specific taxon by NCBI ID or scientific name.",
            operation_id="get_taxon",
            response_description="Returns a `TaxonResponse` object containing detailed information for the specified taxon, such as scientific name, common name, rank, and lineage.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def get_taxon(
            taxon_id: str = FastApiPath(..., description="Taxon identifier (NCBI ID or scientific name)", examples=["9606", "Homo sapiens"]),
            request: Optional[TaxonRequest] = None # Made request body optional
        ):
            """Get taxon information by ID"""
            _fields = request.fields if request and request.fields else "all"
            _email = request.email if request and request.email else None
            _as_dataframe = request.as_dataframe if request else False
            _df_index = request.df_index if request else True            
            
            fields_list = _fields.split(",") if _fields != "all" else None
            
            async with TaxonClientAsync() as client:
                result = await client.gettaxon(
                    taxon_id,
                    fields=fields_list,
                    email=_email,
                    as_dataframe=_as_dataframe,
                    df_index=_df_index
                )
                if result is None:
                    raise HTTPException(status_code=404, detail=f"Taxon '{taxon_id}' not found")
                # Convert to the local TaxonResponse model to ensure correct serialization
                return TaxonResponse.model_validate(result.model_dump(by_alias=True))

        @self.post(
            "/taxons",
            response_model=List[TaxonResponse],
            tags=["taxons"],
            summary="Fetch multiple taxa by a comma-separated list of NCBI IDs or scientific names.",
            operation_id="get_taxons",
            response_description="Returns a list of `TaxonResponse` objects, each containing detailed information (e.g., scientific name, rank) for the corresponding requested taxon ID/name.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def get_taxons(request: TaxonsRequest):
            """Get information for multiple taxa"""
            async with TaxonClientAsync() as client:
                results_from_client = await client.gettaxons(
                    request.taxon_ids.split(","),
                    fields=request.fields.split(",") if request.fields else None,
                    email=request.email,
                    as_dataframe=request.as_dataframe,
                    df_index=request.df_index
                )
                # Convert each item to the local TaxonResponse model
                return [TaxonResponse.model_validate(item.model_dump(by_alias=True)) for item in results_from_client]

        @self.post(
            "/taxon/query",
            response_model=TaxonQueryResponse,
            tags=["taxons"],
            summary="Search taxa via Lucene query (e.g., scientific name, rank), returning taxon details and query metadata.",
            operation_id="query_taxons",
            response_description="Provides a list of `TaxonResponse` objects (e.g., scientific name, rank) as 'hits', along with 'total' count, 'max_score', and 'took' time.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def query_taxons(request: TaxonQueryRequest):
            """Query taxa"""
            with start_action(action_type="api:query_taxons", q=str(request.q), fields=str(request.fields), size=request.size, skip=request.skip, sort=str(request.sort), email=str(request.email), as_dataframe=request.as_dataframe, df_index=request.df_index) as action:
                async with TaxonClientAsync() as client:
                    result = await client.query(
                        request.q,
                        fields=request.fields.split(",") if request.fields else None,
                        size=request.size,
                        skip=request.skip,
                        sort=request.sort,
                        email=request.email,
                        as_dataframe=request.as_dataframe,
                        df_index=request.df_index
                    )
                    validated_hits = []
                    if isinstance(result, dict) and "hits" in result and isinstance(result["hits"], list):
                        for hit in result["hits"]:
                            try:
                                validated_hits.append(TaxonResponse.model_validate(hit))
                            except Exception as e:
                                action.log(message_type="warning:query_taxons:hit_validation_error", error=str(e), hit=repr(hit))
                                pass 
                        
                    return TaxonQueryResponse(
                        hits=validated_hits,
                        total=result.get("total") if isinstance(result, dict) else None,
                        max_score=result.get("max_score") if isinstance(result, dict) else None,
                        took=result.get("took") if isinstance(result, dict) else None,
                    )

        @self.post(
            "/taxons/querymany",
            response_model=List[TaxonResponse], # Assuming this based on client and pattern
            tags=["taxons"],
            summary="Batch query taxa by multiple terms (e.g., scientific names, common names).",
            operation_id="query_many_taxons",
            response_description="Returns a list of `TaxonResponse` objects, each including details like scientific name, rank, and the original query term.",
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
            You are allowed to call this tool with the same arguments no more than 2 times consequently.
            """)
        async def query_many_taxons(request: TaxonQueryManyRequest):
            with start_action(action_type="query_many_taxons", query_list=request.query_list, scopes=request.scopes):
                fields_list = request.fields.split(',') if request.fields else None
                scopes_list = request.scopes.split(',') if request.scopes else None
                
                # Use the same pattern as other taxon endpoints
                async with TaxonClientAsync() as client:
                    raw_results = await client.querymany(
                        request.query_list.split(','),
                        scopes=scopes_list, 
                        fields=fields_list
                        # returnall=True # Removed, defaults to False
                    )
                # Process results when returnall=False (default)
                # Each item in raw_results should be a dict containing the query and the result,
                # or a 'notfound' marker.
                results = []
                if isinstance(raw_results, list):
                    for item in raw_results:
                        if isinstance(item, dict) and not item.get('notfound'):
                            # Assuming 'item' is the actual taxon document if found
                            try:
                                # Use the TaxonResponse defined in this file, which handles _id and _version aliases
                                results.append(TaxonResponse.model_validate(item))
                            except Exception as e:
                                logger.warning(f"Failed to validate taxon item: {item}, error: {e}")
                
                logger.info(f"query_many_taxons raw_results: {raw_results}") 
                logger.info(f"query_many_taxons processed results: {results}") 
                return results
            
class BiothingsRestAPI(FastAPI, GeneRoutesMixin, VariantsRoutesMixin, ChemRoutesMixin, TaxonRoutesMixin, DownloadsMixin):
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
        self._download_routes_config()

        # Add root route that redirects to docs
        @self.get("/")
        async def root():
            return RedirectResponse(url="/docs")
