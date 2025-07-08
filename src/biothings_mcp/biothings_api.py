#!/usr/bin/env python3
"""Biothings API tools for MCP server - converts biological data APIs to MCP tools."""

import logging
from typing import Optional, List, Dict, Any, Union

from pydantic import BaseModel, Field, ConfigDict
from biothings_typed_client.genes import GeneClientAsync, GeneResponse
from biothings_typed_client.variants import VariantClientAsync, VariantResponse
from biothings_typed_client.chem import ChemClientAsync, ChemResponse
from biothings_typed_client.taxons import TaxonClientAsync, TaxonResponse as BaseClientTaxonResponse
from eliot import start_action

from biothings_mcp.download_api import DownloadTools

# Setup logger for this module
logger = logging.getLogger(__name__)

# Custom TaxonResponse model making version field optional
class TaxonResponse(BaseClientTaxonResponse):
    """Response model for taxon information with version field as optional"""
    model_config = ConfigDict(extra='allow')
    
    # Override parent's fields to set new serialization_alias and modify type (for version)
    # Also ensure validation_alias is present for input mapping
    id: str = Field(validation_alias='_id', serialization_alias='_id', description="Taxon identifier")
    version: Optional[int] = Field(default=1, validation_alias='_version', serialization_alias='_version', description="Version number")

# Request Body Models for MCP tools
class GeneQueryRequest(BaseModel):
    q: str = Field(..., description="Query string following Lucene syntax.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
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
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    species: Optional[str] = Field(None, description="Filter by species.")
    email: Optional[str] = Field(None, description="User email for tracking.")
    as_dataframe: bool = Field(False, description="Return as DataFrame.")
    df_index: bool = Field(True, description="Index DataFrame by matched query term.")
    size: int = Field(10, description="Maximum number of hits per query term.")

class GeneRequest(BaseModel):
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    species: Optional[str] = Field(None, description="Specify species.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class GenesRequest(BaseModel):
    gene_ids: str = Field(..., description="Comma-separated list of gene IDs.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    species: Optional[str] = Field(None, description="Filter by species.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class VariantQueryRequest(BaseModel):
    q: str = Field(..., description="Query string following Lucene syntax.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    size: int = Field(10, description="Maximum number of hits.")
    skip: int = Field(0, description="Number of hits to skip.")
    sort: Optional[str] = Field(None, description="Comma-separated fields to sort on.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class VariantQueryManyRequest(BaseModel):
    query_list: str = Field(..., description="Comma-separated list of query terms.")
    scopes: Optional[str] = Field(None, description="Comma-separated list of fields to search against.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class VariantRequest(BaseModel):
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class VariantsRequest(BaseModel):
    variant_ids: str = Field(..., description="Comma-separated list of variant IDs.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class ChemQueryRequest(BaseModel):
    q: str = Field(..., description="Query string.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    size: int = Field(10, description="Maximum number of results.")
    skip: int = Field(0, description="Number of results to skip.")
    sort: Optional[str] = Field(None, description="Sort field.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class ChemQueryManyRequest(BaseModel):
    query_list: str = Field(..., description="Comma-separated list of query strings.")
    scopes: Optional[str] = Field(None, description="Comma-separated list of fields to search in.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class ChemRequest(BaseModel):
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class ChemsRequest(BaseModel):
    chem_ids: str = Field(..., description="Comma-separated list of chemical IDs.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class TaxonRequest(BaseModel):
    fields: str = Field("all", description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class TaxonsRequest(BaseModel):
    taxon_ids: str = Field(..., description="Comma-separated list of taxon IDs.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class TaxonQueryRequest(BaseModel):
    q: str = Field(..., description="Query string.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    size: int = Field(10, description="Maximum number of results.")
    skip: int = Field(0, description="Number of results to skip.")
    sort: Optional[str] = Field(None, description="Sort field.")
    email: Optional[str] = Field(None, description="User email for tracking.")

class TaxonQueryManyRequest(BaseModel):
    query_list: str = Field(..., description="Comma-separated list of query strings.")
    scopes: Optional[str] = Field(None, description="Comma-separated list of fields to search in.")
    fields: Optional[str] = Field("all", description="Comma-separated list of fields to return.")
    email: Optional[str] = Field(None, description="User email for tracking.")

# Response models
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

class ChemQueryResponse(BaseModel):
    hits: List[ChemResponse]
    total: Optional[int] = None
    max_score: Optional[float] = None
    took: Optional[int] = None

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

# Tool Handlers for MCP
class GeneTools:
    """Handler for gene-related MCP tools."""
    
    def __init__(self, mcp_server, prefix: str = ""):
        self.mcp_server = mcp_server
        self.prefix = prefix
        self.gene_client = GeneClientAsync()
    
    async def query_genes(self, q: str, fields: Optional[str] = "all", size: int = 10, skip: int = 0, 
                   sort: Optional[str] = None, species: Optional[str] = None, 
                   email: Optional[str] = None) -> QueryResponse:
        """Search genes via Lucene query, returning gene details and query metadata.
        
        **IMPORTANT:** This endpoint requires structured queries using specific field names. 
        Simple natural language queries like "CDK2 gene" or "human kinase" will **NOT** work.
        You **MUST** specify the field you are querying, e.g., `symbol:CDK2`, `name:"cyclin-dependent kinase 2"`, `taxid:9606`.
        Use this tool when you need to *search* for genes based on criteria, not when you already know the specific gene ID.
        If you know the exact Entrez or Ensembl ID, use the `get_gene` tool instead for faster retrieval.
        
        **Supported Query Features (based on Lucene syntax):**
        1. Simple Term Queries: `q=cdk2` (Searches across default fields)
        2. Fielded Queries: `q=symbol:CDK2`, `q=name:"cyclin-dependent kinase 2"`
        3. Range Queries: `q=taxid:[9606 TO 10090]`
        4. Boolean Queries: `q=symbol:CDK2 AND taxid:9606`
        5. Wildcard Queries: `q=symbol:CDK*`
        
        Returns gene details including symbol, name, taxid, and entrezgene.
        """
        with start_action(action_type="query_genes", q=q, size=size) as action:
            async with self.gene_client as client:
                raw_result = await client.query(
                    q=q, fields=fields, size=size, skip=skip, sort=sort, 
                    species=species, email=email
                )
            # Convert raw dict to QueryResponse model
            if isinstance(raw_result, dict):
                result = QueryResponse(
                    hits=[GeneResponse.model_validate(hit) for hit in raw_result.get('hits', [])],
                    total=raw_result.get('total'),
                    max_score=raw_result.get('max_score'), 
                    took=raw_result.get('took')
                )
            else:
                result = raw_result
            action.add_success_fields(hits_count=len(result.hits))
            return result
    
    async def query_many_genes(self, query_list: str, scopes: Optional[str] = "entrezgene,ensemblgene,retired",
                       fields: Optional[str] = "all", species: Optional[str] = None,
                       email: Optional[str] = None, size: int = 10) -> List[GeneResponse]:
        """Batch query genes by multiple terms, returning a list of gene details.
        
        Perform multiple gene searches in a single request using a comma-separated list of query terms.
        Unlike `query_genes`, the `query_list` parameter takes multiple **terms** (like gene IDs, symbols, names) 
        rather than full query strings. The `scopes` parameter defines which fields these terms should be searched against.
        
        **Endpoint Usage:**
        - Query multiple symbols: `query_list=CDK2,BRCA1` with `scopes=symbol`
        - Query multiple Entrez IDs: `query_list=1017,672` with `scopes=entrezgene`
        - Query mixed IDs/symbols: `query_list=CDK2,672` with `scopes=symbol,entrezgene`
        """
        with start_action(action_type="query_many_genes", query_list=query_list) as action:
            async with self.gene_client as client:
                raw_result = await client.querymany(
                    query_list=query_list.split(','), scopes=scopes, fields=fields,
                    species=species, email=email, size=size
                )
            # Convert raw dicts to GeneResponse models, handling 'notfound' entries
            if isinstance(raw_result, list):
                result = []
                for item in raw_result:
                    if isinstance(item, dict) and not item.get('notfound', False):
                        try:
                            result.append(GeneResponse.model_validate(item))
                        except Exception:
                            # Skip invalid entries
                            continue
            else:
                result = raw_result
            action.add_success_fields(results_count=len(result))
            return result
    
    async def get_gene(self, gene_id: str, fields: Optional[str] = "all", species: Optional[str] = None,
                email: Optional[str] = None) -> GeneResponse:
        """Fetch a specific gene by Entrez or Ensembl ID.
        
        Retrieves detailed information for a **single, specific gene** using its exact known identifier.
        **This is the preferred tool over `query_genes` for fetching a specific gene when you already know 
        its standard ID (Entrez or Ensembl) and don't need complex search filters.**
        
        **Supported Identifiers:**
        - Entrez Gene ID: e.g., `1017`
        - Ensembl Gene ID: e.g., `ENSG00000123374`
        """
        with start_action(action_type="get_gene", gene_id=gene_id) as action:
            async with self.gene_client as client:
                result = await client.getgene(
                    gene_id=gene_id, fields=fields, species=species, email=email
                )
            action.add_success_fields(gene_found=True)
            return result
    
    async def get_genes(self, gene_ids: str, fields: Optional[str] = "all", species: Optional[str] = None,
                 email: Optional[str] = None) -> List[GeneResponse]:
        """Fetch multiple genes by a comma-separated list of Entrez or Ensembl IDs.
        
        Retrieves detailed information for **multiple specific genes** in a single request using their exact known identifiers.
        **This is the preferred tool over `query_many_genes` for fetching multiple specific genes when you already know 
        their standard IDs (Entrez, Ensembl) and don't need complex search filters.**
        
        **Input Format:** Accepts a comma-separated list of gene IDs (Entrez or Ensembl).
        **Example:** `gene_ids=1017,1018` or `gene_ids=ENSG00000123374,ENSG00000134057`
        """
        with start_action(action_type="get_genes", gene_ids=gene_ids) as action:
            async with self.gene_client as client:
                result = await client.getgenes(
                    gene_ids=gene_ids.split(','), fields=fields, species=species, email=email
                )
            action.add_success_fields(genes_count=len(result))
            return result
    
    async def get_gene_metadata(self, email: Optional[str] = None) -> MetadataResponse:
        """Retrieve MyGene.info database metadata including stats and fields.
        
        Retrieve metadata about the underlying MyGene.info gene annotation database, **NOT** information about specific genes.
        **Use this tool ONLY to understand the database itself** (e.g., to discover available fields, check data versions, 
        or get overall statistics). It **CANNOT** be used to find or retrieve data for any particular gene.
        
        **Returned Information:**
        - `stats`: Database statistics (e.g., total number of genes)
        - `fields`: Available gene annotation fields and their data types
        - `index`: Information about the backend data index
        - `version`: Data version information
        """
        with start_action(action_type="get_gene_metadata") as action:
            async with self.gene_client as client:
                raw_result = await client.metadata()
            # Convert raw dict to MetadataResponse model
            if isinstance(raw_result, dict):
                result = MetadataResponse(
                    stats=raw_result.get('stats', {}),
                    fields=raw_result.get('fields'),
                    index=raw_result.get('index'),
                    version=raw_result.get('version')
                )
            else:
                result = raw_result
            action.add_success_fields(metadata_retrieved=True)
            return result
    
    def register_tools(self):
        """Register gene-related MCP tools."""
        
        self.mcp_server.tool(
            name=f"{self.prefix}query_genes",
            description=self.query_genes.__doc__
        )(self.query_genes)
        
        self.mcp_server.tool(
            name=f"{self.prefix}query_many_genes",
            description=self.query_many_genes.__doc__
        )(self.query_many_genes)
        
        self.mcp_server.tool(
            name=f"{self.prefix}get_gene",
            description=self.get_gene.__doc__
        )(self.get_gene)
        
        self.mcp_server.tool(
            name=f"{self.prefix}get_genes",
            description=self.get_genes.__doc__
        )(self.get_genes)
        
        self.mcp_server.tool(
            name=f"{self.prefix}get_gene_metadata",
            description=self.get_gene_metadata.__doc__
        )(self.get_gene_metadata)

class VariantTools:
    """Handler for variant-related MCP tools."""
    
    def __init__(self, mcp_server, prefix: str = ""):
        self.mcp_server = mcp_server
        self.prefix = prefix
        self.variant_client = VariantClientAsync()
    
    async def query_variants(self, q: str, fields: Optional[str] = "all", size: int = 10, skip: int = 0,
                      sort: Optional[str] = None, email: Optional[str] = None) -> VariantQueryResponse:
        """Search variants via Lucene query (e.g., rsID, gene name), returning variant details and query metadata.
        
        Search for variants using a query string with various filtering options, leveraging the MyVariant.info API.
        **Use this tool for *searching* variants based on criteria.** 
        If you already know the exact variant ID (HGVS, rsID), use the `get_variant` tool for faster direct retrieval.

        **Supported Query Features (Lucene syntax):**
        1. Simple Queries: `q=rs58991260` (Find by rsID)
        2. Fielded Queries: `q=dbsnp.vartype:snp`, `q=dbnsfp.polyphen2.hdiv.pred:(D P)`
        3. Range Queries: `q=dbnsfp.polyphen2.hdiv.score:>0.99`
        4. Wildcard Queries: `q=dbnsfp.genename:CDK*`
        5. Boolean Queries: `q=_exists_:dbsnp AND dbsnp.vartype:snp`
        6. Genomic Interval Queries: `q=chr1:69000-70000`
        """
        with start_action(action_type="query_variants", q=q, size=size) as action:
            async with self.variant_client as client:
                raw_result = await client.query(
                    q=q, fields=fields, size=size, skip=skip, sort=sort, email=email
                )
            # Convert raw dict to VariantQueryResponse model
            if isinstance(raw_result, dict):
                result = VariantQueryResponse(
                    hits=[VariantResponse.model_validate(hit) for hit in raw_result.get('hits', [])],
                    total=raw_result.get('total'),
                    max_score=raw_result.get('max_score'),
                    took=raw_result.get('took')
                )
            else:
                result = raw_result
            action.add_success_fields(hits_count=len(result.hits))
            return result
    
    async def query_many_variants(self, query_list: str, scopes: Optional[str] = None, fields: Optional[str] = "all",
                           email: Optional[str] = None) -> List[VariantResponse]:
        """Batch query variants by multiple identifiers (e.g., rsIDs, HGVS IDs).
        
        Perform multiple variant queries in a single request using a comma-separated list of variant identifiers.
        This tool takes multiple **terms** (like rsIDs, HGVS IDs) in `query_list` and searches for them within the specified `scopes`.
        
        **Endpoint Usage:**
        - Query multiple rsIDs: `query_list=rs58991260,rs2500` with `scopes=dbsnp.rsid`
        - Query multiple HGVS IDs: `query_list=chr7:g.140453134T>C,chr1:g.69511A>G`
        - Query mixed IDs: `query_list=rs58991260,chr1:g.69511A>G` with `scopes=dbsnp.rsid,_id`
        """
        with start_action(action_type="query_many_variants", query_list=query_list) as action:
            async with self.variant_client as client:
                raw_result = await client.querymany(
                    query_list=query_list.split(','), scopes=scopes, fields=fields, email=email
                )
            # Convert raw dicts to VariantResponse models, handling 'notfound' entries
            if isinstance(raw_result, list):
                result = []
                for item in raw_result:
                    if isinstance(item, dict) and not item.get('notfound', False):
                        try:
                            result.append(VariantResponse.model_validate(item))
                        except Exception:
                            # Skip invalid entries
                            continue
            else:
                result = raw_result
            action.add_success_fields(results_count=len(result))
            return result
    
    async def get_variant(self, variant_id: str, fields: Optional[str] = "all", email: Optional[str] = None) -> VariantResponse:
        """Fetch a specific variant by HGVS or rsID.
        
        Retrieves detailed annotation data for a **single, specific variant** using its identifier.
        **This is the preferred tool over `query_variants` for fetching a specific variant when you already know 
        its standard ID (HGVS or rsID) and don't need complex search filters.**
        
        **Supported Identifiers:**
        - HGVS ID (e.g., `chr7:g.140453134T>C`). *Note: MyVariant.info primarily uses hg19-based HGVS IDs.*
        - dbSNP rsID (e.g., `rs58991260`)
        """
        with start_action(action_type="get_variant", variant_id=variant_id) as action:
            async with self.variant_client as client:
                # When requesting specific fields, we need to ensure required fields are included
                if fields:
                    # Add required fields for VariantResponse validation
                    fields_list = fields.split(',') if isinstance(fields, str) else fields
                    required_fields = ['_id', 'chrom', 'vcf.alt', 'vcf.position', 'vcf.ref']
                    for req_field in required_fields:
                        if req_field not in fields_list:
                            fields_list.append(req_field)
                    fields = ','.join(fields_list)
                
                result = await client.getvariant(
                    variant_id=variant_id, fields=fields, email=email
                )
            action.add_success_fields(variant_found=True)
            return result
    
    async def get_variants(self, variant_ids: str, fields: Optional[str] = "all", email: Optional[str] = None) -> List[VariantResponse]:
        """Fetch multiple variants by a comma-separated list of HGVS or rsIDs.
        
        Retrieves annotation data for **multiple specific variants** in a single request using their identifiers.
        **This is the preferred tool over `query_many_variants` for fetching multiple specific variants when you already know 
        their standard IDs (HGVS or rsID).**

        **Input Format:** Accepts a comma-separated list of variant IDs (HGVS or dbSNP rsIDs).
        **Examples:** `variant_ids=chr7:g.140453134T>C,chr1:g.69511A>G` or `variant_ids=rs58991260,rs2500`
        """
        with start_action(action_type="get_variants", variant_ids=variant_ids) as action:
            async with self.variant_client as client:
                raw_result = await client.getvariants(
                    variant_ids=variant_ids.split(','), fields=fields, email=email
                )
            # Filter out 'notfound' entries and convert to VariantResponse models
            result = []
            if isinstance(raw_result, list):
                for item in raw_result:
                    if isinstance(item, dict) and not item.get('notfound', False):
                        try:
                            result.append(VariantResponse.model_validate(item))
                        except Exception:
                            # Skip invalid entries
                            continue
                    elif hasattr(item, 'model_validate'):
                        # Already a VariantResponse object
                        result.append(item)
            else:
                result = raw_result
            action.add_success_fields(variants_count=len(result))
            return result
    
    def register_tools(self):
        """Register variant-related MCP tools."""
        
        self.mcp_server.tool(
            name=f"{self.prefix}query_variants",
            description=self.query_variants.__doc__
        )(self.query_variants)
        
        self.mcp_server.tool(
            name=f"{self.prefix}query_many_variants",
            description=self.query_many_variants.__doc__
        )(self.query_many_variants)
        
        self.mcp_server.tool(
            name=f"{self.prefix}get_variant",
            description=self.get_variant.__doc__
        )(self.get_variant)
        
        self.mcp_server.tool(
            name=f"{self.prefix}get_variants",
            description=self.get_variants.__doc__
        )(self.get_variants)

class ChemTools:
    """Handler for chemical-related MCP tools."""
    
    def __init__(self, mcp_server, prefix: str = ""):
        self.mcp_server = mcp_server
        self.prefix = prefix
        self.chem_client = ChemClientAsync()
    
    async def query_chems(self, q: str, fields: Optional[str] = "all", size: int = 10, skip: int = 0,
                   sort: Optional[str] = None, email: Optional[str] = None) -> ChemQueryResponse:
        """Search chemical compounds via Lucene query (e.g., name, formula), returning compound details and query metadata.
        
        Search for chemical compounds using a query string with various filtering options.
        
        **Supported Query Features:**
        1. Simple Queries: "C6H12O6" - Find compounds with molecular formula, "glucose" - Find compounds with name
        2. Fielded Queries: "pubchem.molecular_formula:C6H12O6", "pubchem.molecular_weight:[100 TO 200]"
        3. Range Queries: "pubchem.xlogp:>2", "pubchem.topological_polar_surface_area:[50 TO 100]"
        4. Boolean Queries: "pubchem.hydrogen_bond_donor_count:>2 AND pubchem.hydrogen_bond_acceptor_count:>4"
        
        Returns compound details including PubChem data like formula, weight, and XLogP.
        """
        with start_action(action_type="query_chems", q=q, size=size) as action:
            async with self.chem_client as client:
                raw_result = await client.query(
                    q=q, fields=fields, size=size, skip=skip, sort=sort, email=email
                )
            # Convert raw dict to ChemQueryResponse model
            if isinstance(raw_result, dict):
                result = ChemQueryResponse(
                    hits=[ChemResponse.model_validate(hit) for hit in raw_result.get('hits', [])],
                    total=raw_result.get('total'),
                    max_score=raw_result.get('max_score'),
                    took=raw_result.get('took')
                )
            else:
                result = raw_result
            action.add_success_fields(hits_count=len(result.hits))
            return result
    
    async def query_many_chems(self, query_list: str, scopes: Optional[str] = None, fields: Optional[str] = "all",
                        email: Optional[str] = None) -> List[ChemResponse]:
        """Batch query chemical compounds by multiple terms (e.g., names, InChIKeys).
        
        Perform multiple chemical queries in a single request.
        
        **Supported Usage:**
        1. Multiple Query Types: ["C6H12O6", "C12H22O11"] (formulas), ["glucose", "sucrose"] (names)
        2. Field Scoping: Search in specific fields using scopes parameter
        3. Result Filtering: Return specific fields using fields parameter
        """
        with start_action(action_type="query_many_chems", query_list=query_list) as action:
            async with self.chem_client as client:
                raw_result = await client.querymany(
                    query_list=query_list.split(','), scopes=scopes, fields=fields, email=email
                )
            # Convert raw dicts to ChemResponse models, handling 'notfound' entries
            if isinstance(raw_result, list):
                result = []
                for item in raw_result:
                    if isinstance(item, dict) and not item.get('notfound', False):
                        try:
                            result.append(ChemResponse.model_validate(item))
                        except Exception:
                            # Skip invalid entries
                            continue
            else:
                result = raw_result
            action.add_success_fields(results_count=len(result))
            return result
    
    async def get_chem(self, chem_id: str, fields: Optional[str] = "all", email: Optional[str] = None) -> ChemResponse:
        """Fetch a specific chemical compound by ID (e.g., InChIKey, PubChem CID).
        
        Retrieves detailed information about a specific chemical compound using its identifier.
        
        **Supported ID formats:**
        - InChIKey: "KTUFNOKKBVMGRW-UHFFFAOYSA-N" (Glucose)
        - PubChem CID: "5793" (Glucose)
        - SMILES: "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O"
        
        Returns comprehensive chemical information including structural properties, physical properties, 
        chemical properties, stereochemistry information, and chemical identifiers.
        """
        with start_action(action_type="get_chem", chem_id=chem_id) as action:
            async with self.chem_client as client:
                result = await client.getchem(
                    chem_id=chem_id, fields=fields, email=email
                )
            action.add_success_fields(chem_found=True)
            return result
    
    async def get_chems(self, chem_ids: str, fields: Optional[str] = "all", email: Optional[str] = None) -> List[ChemResponse]:
        """Fetch multiple chemical compounds by a comma-separated list of IDs.
        
        Retrieves information for multiple chemical compounds in a single request.
        
        **Input Format:** Accepts comma-separated list of chemical IDs in various formats:
        - InChIKeys: "KTUFNOKKBVMGRW-UHFFFAOYSA-N,XEFQLINVKFYRCS-UHFFFAOYSA-N"
        - PubChem CIDs: "5793,5281"
        - Mixed formats: "KTUFNOKKBVMGRW-UHFFFAOYSA-N,5281"
        """
        with start_action(action_type="get_chems", chem_ids=chem_ids) as action:
            async with self.chem_client as client:
                result = await client.getchems(
                    chem_ids=chem_ids.split(','), fields=fields, email=email
                )
            action.add_success_fields(chems_count=len(result))
            return result
    
    def register_tools(self):
        """Register chemical-related MCP tools."""
        
        self.mcp_server.tool(
            name=f"{self.prefix}query_chems",
            description=self.query_chems.__doc__
        )(self.query_chems)
        
        self.mcp_server.tool(
            name=f"{self.prefix}query_many_chems",
            description=self.query_many_chems.__doc__
        )(self.query_many_chems)
        
        self.mcp_server.tool(
            name=f"{self.prefix}get_chem",
            description=self.get_chem.__doc__
        )(self.get_chem)
        
        self.mcp_server.tool(
            name=f"{self.prefix}get_chems",
            description=self.get_chems.__doc__
        )(self.get_chems)

class TaxonTools:
    """Handler for taxon-related MCP tools."""
    
    def __init__(self, mcp_server, prefix: str = ""):
        self.mcp_server = mcp_server
        self.prefix = prefix
        self.taxon_client = TaxonClientAsync()
    
    async def get_taxon(self, taxon_id: str, fields: str = "all", email: Optional[str] = None) -> TaxonResponse:
        """Fetch a specific taxon by NCBI ID or scientific name.
        
        Retrieves detailed information about a specific taxon using its identifier.
        
        **Supported Identifiers:**
        - NCBI ID: 9606 (Homo sapiens)
        - Scientific name: "Homo sapiens"
        
        Returns comprehensive taxon information including basic information (ID, scientific name, common name),
        taxonomic classification (rank, parent taxon), lineage information, alternative names and authorities,
        and gene data availability.
        """
        with start_action(action_type="get_taxon", taxon_id=taxon_id) as action:
            async with self.taxon_client as client:
                result = await client.gettaxon(
                    taxon_id=taxon_id, fields=fields, email=email
                )
            action.add_success_fields(taxon_found=True)
            return TaxonResponse.model_validate(result.model_dump(by_alias=True))
    
    async def get_taxons(self, taxon_ids: str, fields: Optional[str] = "all", email: Optional[str] = None) -> List[TaxonResponse]:
        """Fetch multiple taxa by a comma-separated list of NCBI IDs or scientific names.
        
        Retrieves information for multiple taxa in a single request.
        
        **Input Format:** Accepts comma-separated list of taxon IDs (either NCBI IDs or scientific names).
        **Examples:**
        - Multiple NCBI IDs: "9606,10090" (Homo sapiens and Mus musculus)
        - Multiple scientific names: "Homo sapiens,Mus musculus"
        - Mixed IDs: "9606,Mus musculus" (Homo sapiens by NCBI ID and Mus musculus by name)
        """
        with start_action(action_type="get_taxons", taxon_ids=taxon_ids) as action:
            async with self.taxon_client as client:
                # Get individual taxons since gettaxons doesn't exist
                taxon_id_list = taxon_ids.split(',')
                results_from_client = []
                for taxon_id in taxon_id_list:
                    try:
                        taxon_result = await client.gettaxon(
                            taxon_id=taxon_id.strip(), fields=fields, email=email
                        )
                        results_from_client.append(taxon_result)
                    except Exception:
                        # Skip failed taxons
                        continue
            # Convert each item to the local TaxonResponse model
            result = [TaxonResponse.model_validate(item.model_dump(by_alias=True)) for item in results_from_client]
            action.add_success_fields(taxons_count=len(result))
            return result
    
    async def query_taxons(self, q: str, fields: Optional[str] = "all", size: int = 10, skip: int = 0,
                    sort: Optional[str] = None, email: Optional[str] = None) -> TaxonQueryResponse:
        """Search taxa via Lucene query (e.g., scientific name, rank), returning taxon details and query metadata.
        
        Search for taxa using a query string with various filtering options.
        
        **Supported Query Features:**
        1. Simple Queries: "scientific_name:Homo sapiens", "common_name:human"
        2. Fielded Queries: "rank:species", "parent_taxid:9606", "has_gene:true"
        3. Range Queries: "taxid:[9606 TO 10090]", "lineage:>9606"
        4. Boolean Queries: "rank:species AND has_gene:true", "scientific_name:Homo* AND NOT rank:genus"
        5. Wildcard Queries: "scientific_name:Homo*", "common_name:*mouse*"
        """
        with start_action(action_type="query_taxons", q=q, size=size) as action:
            async with self.taxon_client as client:
                result = await client.query(
                    q=q, fields=fields, size=size, skip=skip, sort=sort, email=email
                )
            validated_hits = []
            if isinstance(result, dict) and "hits" in result and isinstance(result["hits"], list):
                for hit in result["hits"]:
                    try:
                        validated_hits.append(TaxonResponse.model_validate(hit))
                    except Exception as e:
                        logger.warning(f"Failed to validate taxon hit: {hit}, error: {e}")
                        pass 
            
            query_response = TaxonQueryResponse(
                hits=validated_hits,
                total=result.get("total") if isinstance(result, dict) else None,
                max_score=result.get("max_score") if isinstance(result, dict) else None,
                took=result.get("took") if isinstance(result, dict) else None,
            )
            action.add_success_fields(hits_count=len(validated_hits))
            return query_response
    
    async def query_many_taxons(self, query_list: str, scopes: Optional[str] = None, fields: Optional[str] = "all",
                         email: Optional[str] = None) -> List[TaxonResponse]:
        """Batch query taxa by multiple terms (e.g., scientific names, common names).
        
        Perform multiple taxon queries in a single request.
        
        **Supported Usage:**
        1. Multiple Query Types: ["Homo sapiens", "Mus musculus"] (scientific names), ["human", "mouse"] (common names)
        2. Field Scoping: Search in specific fields using scopes parameter: ["scientific_name", "common_name"]
        3. Result Filtering: Return specific fields using fields parameter: ["scientific_name", "common_name", "rank"]
        """
        with start_action(action_type="query_many_taxons", query_list=query_list) as action:
            fields_list = fields.split(',') if fields else None
            scopes_list = scopes.split(',') if scopes else None
            
            async with self.taxon_client as client:
                raw_results = await client.querymany(
                    query_list.split(','),
                    scopes=scopes_list, 
                    fields=fields_list
                )
            
            # Process results when returnall=False (default)
            results = []
            if isinstance(raw_results, list):
                for item in raw_results:
                    if isinstance(item, dict) and not item.get('notfound'):
                        try:
                            results.append(TaxonResponse.model_validate(item))
                        except Exception as e:
                            logger.warning(f"Failed to validate taxon item: {item}, error: {e}")
            
            action.add_success_fields(results_count=len(results))
            return results
    
    def register_tools(self):
        """Register taxon-related MCP tools."""
        
        self.mcp_server.tool(
            name=f"{self.prefix}get_taxon",
            description=self.get_taxon.__doc__
        )(self.get_taxon)
        
        self.mcp_server.tool(
            name=f"{self.prefix}get_taxons",
            description=self.get_taxons.__doc__
        )(self.get_taxons)
        
        self.mcp_server.tool(
            name=f"{self.prefix}query_taxons",
            description=self.query_taxons.__doc__
        )(self.query_taxons)
        
        self.mcp_server.tool(
            name=f"{self.prefix}query_many_taxons",
            description=self.query_many_taxons.__doc__
        )(self.query_many_taxons)


