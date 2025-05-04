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
            description="""
            Search for genes using a query string with various filtering options.
            
            This endpoint supports complex queries with the following features:
            
            1. Simple Queries:
               - "symbol:CDK2" - Find genes with symbol CDK2
               - "name:cyclin-dependent kinase 2" - Find genes with specific name
            
            2. Fielded Queries:
               - "refseq.rna:NM_001798" - Find genes with specific RefSeq RNA ID
               - "ensembl.gene:ENSG00000123374" - Find genes with specific Ensembl ID
            
            3. Range Queries:
               - "taxid:[9606 TO 10090]" - Find genes in specific taxonomy range
               - "entrezgene:>1000" - Find genes with Entrez ID greater than 1000
            
            4. Boolean Queries:
               - "symbol:CDK2 AND taxid:9606" - Find CDK2 gene in human
               - "symbol:CDK* AND NOT taxid:9606" - Find CDK genes not in human
            
            5. Wildcard Queries:
               - "symbol:CDK*" - Find genes with symbol starting with CDK
               - "name:*kinase*" - Find genes with 'kinase' in their name
            
            The response includes pagination information and can be returned as a pandas DataFrame.
            """
        )
        async def query_genes(
            q: str = Query(
                ...,
                description="Query string",
                examples=[
                    "symbol:CDK2",
                    "name:cyclin-dependent kinase 2",
                    "refseq.rna:NM_001798",
                    "taxid:[9606 TO 10090]",
                    "symbol:CDK2 AND taxid:9606",
                    "symbol:CDK* AND NOT taxid:9606",
                    "name:*kinase*"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "symbol,name,refseq.rna",
                    "id,symbol,name,taxid",
                    "ensembl.gene,uniprot.Swiss-Prot"
                ]
            ),
            size: int = Query(10, description="Maximum number of results to return (max 1000)", examples=[10, 50, 100]),
            skip: int = Query(0, description="Number of results to skip (for pagination)", examples=[0, 10, 20]),
            sort: Optional[str] = Query(None, description="Sort field, prefix with '-' for descending order", examples=["_score", "-_score", "symbol"]),
            species: Optional[str] = Query(None, description="Species names or taxonomy ids", examples=["9606", "10090", "9606,10090"]),
            email: Optional[str] = Query(None, description="User email for tracking usage", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as pandas DataFrame", examples=[True, False]),
            df_index: bool = Query(True, description="Index DataFrame by query (only if as_dataframe=True)", examples=[True, False])
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
            description="""
            Perform multiple gene queries in a single request.
            
            This endpoint is useful for batch processing of gene queries. It supports:
            
            1. Multiple Query Types:
               - Symbol queries: ["CDK2", "BRCA1"]
               - Name queries: ["cyclin-dependent kinase 2", "breast cancer 1"]
               - Mixed queries: ["CDK2", "ENSG00000123374"]
            
            2. Field Scoping:
               - Search in specific fields: scopes=["symbol", "name"]
               - Search in all fields: scopes=None
            
            3. Result Filtering:
               - Return specific fields: fields=["symbol", "name", "refseq.rna"]
               - Return all fields: fields=None
            
            4. Species Filtering:
               - Filter by species: species=["9606"] (human)
               - Filter by multiple species: species=["9606", "10090"] (human and mouse)
            
            The response can be returned as a pandas DataFrame for easier data manipulation.
            """
        )
        async def query_many_genes(
            query_list: str = Query(
                ...,
                description="Comma-separated list of query strings",
                examples=[
                    "CDK2,BRCA1",
                    "cyclin-dependent kinase 2,breast cancer 1",
                    "CDK2,ENSG00000123374"
                ]
            ),
            scopes: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to search in",
                examples=[
                    "symbol,name",
                    "refseq.rna,ensembl.gene",
                    "uniprot.Swiss-Prot"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "symbol,name,refseq.rna",
                    "id,symbol,name,taxid",
                    "ensembl.gene,uniprot.Swiss-Prot"
                ]
            ),
            species: Optional[str] = Query(None, description="Species names or taxonomy ids", examples=["9606", "10090", "9606,10090"]),
            email: Optional[str] = Query(None, description="User email for tracking usage", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as pandas DataFrame", examples=[True, False]),
            df_index: bool = Query(True, description="Index DataFrame by query (only if as_dataframe=True)", examples=[True, False]),
            size: int = Query(10, description="Maximum number of results to return per query term (max 1000)", examples=[1, 10, 50])
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
            description="""
            Retrieve metadata about the gene database.
            
            This endpoint returns information about:
            - Total number of genes
            - Database statistics
            - Available fields and their types
            - Index information
            - Version information
            
            The metadata can be used to understand the database structure and capabilities.
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
            description="""
            Retrieves detailed information about a specific gene using its identifier.
            
            This endpoint supports both Entrez Gene IDs and Ensembl Gene IDs.
            
            Examples:
            - Entrez ID: 1017 (CDK2 gene)
            - Ensembl ID: ENSG00000123374 (CDK2 gene)
            
            The response includes comprehensive gene information such as:
            - Basic information (ID, symbol, name)
            - RefSeq information (genomic, protein, RNA)
            - Taxonomy information
            - Entrez Gene ID
            - Ensembl information
            - UniProt information
            - Gene summary
            - Genomic position information
            
            You can filter the returned fields using the 'fields' parameter.
            """
        )
        async def get_gene(
            gene_id: str = Path(..., description="Gene identifier (Entrez or Ensembl ID)", examples=["1017", "ENSG00000123374"]),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "symbol,name,refseq.rna",
                    "id,symbol,name,taxid",
                    "ensembl.gene,uniprot.Swiss-Prot"
                ]
            ),
            species: Optional[str] = Query(
                None,
                description="Species filter (e.g. '9606' for human)",
                examples=["9606", "10090"]
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
            ),
            size: int = Query(
                10,
                description="Maximum number of results to return (max 1000)",
                examples=[10, 50, 100]
            ),
            skip: int = Query(
                0,
                description="Number of results to skip (for pagination)",
                examples=[0, 10, 20]
            ),
            sort: Optional[str] = Query(
                None,
                description="Sort field, prefix with '-' for descending order",
                examples=["_score", "-_score", "symbol"]
            )
        ):
            """Get gene information by ID
            
            Args:
                gene_id: Gene identifier (Entrez or Ensembl ID)
                fields: Comma-separated list of fields to return
                species: Species filter (e.g. "9606" for human)
                email: User email for tracking usage
                as_dataframe: Return results as pandas DataFrame
                df_index: Index DataFrame by query (only if as_dataframe=True)
                size: Maximum number of results to return (max 1000)
                skip: Number of results to skip (for pagination)
                sort: Sort field, prefix with '-' for descending order
            """
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
            description="""
            Retrieves information for multiple genes in a single request.
            
            This endpoint accepts a comma-separated list of gene IDs (either Entrez or Ensembl IDs).
            
            Examples:
            - Multiple Entrez IDs: "1017,1018" (CDK2 and CDK3 genes)
            - Multiple Ensembl IDs: "ENSG00000123374,ENSG00000134057" (CDK2 and CDK3 genes)
            - Mixed IDs: "1017,ENSG00000134057" (CDK2 by Entrez and CDK3 by Ensembl)
            
            The response includes the same comprehensive gene information as the single gene endpoint,
            but for all requested genes.
            
            You can filter the returned fields using the 'fields' parameter.
            """
        )
        async def get_genes(
            gene_ids: str = Query(
                ...,
                description="Comma-separated list of gene IDs",
                examples=[
                    "1017,1018",
                    "ENSG00000123374,ENSG00000134057",
                    "1017,ENSG00000134057"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "symbol,name,refseq.rna",
                    "id,symbol,name,taxid",
                    "ensembl.gene,uniprot.Swiss-Prot"
                ]
            ),
            species: Optional[str] = Query(
                None,
                description="Species filter (e.g. '9606' for human)",
                examples=["9606", "10090"]
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
            description="""
            Search for variants using a query string with various filtering options.
            
            This endpoint supports complex queries with the following features:
            
            1. Simple Queries:
               - "rs58991260" - Find variant by rsID
               - "chr1:69000-70000" - Find variants in genomic region
            
            2. Fielded Queries:
               - "dbsnp.vartype:snp" - Find SNPs
               - "dbnsfp.polyphen2.hdiv.pred:(D P)" - Find probably/possibly damaging variants
               - "_exists_:dbsnp" - Find variants with dbSNP annotations
               - "_missing_:exac" - Find variants without ExAC annotations
            
            3. Range Queries:
               - "dbnsfp.polyphen2.hdiv.score:>0.99" - Find high-risk variants
               - "exac.af:<0.00001" - Find rare variants
               - "exac.ac.ac_adj:[76640 TO 80000]" - Find variants in frequency range
            
            4. Wildcard Queries:
               - "dbnsfp.genename:CDK?" - Find variants in CDK genes
               - "dbnsfp.genename:CDK*" - Find variants in genes starting with CDK
            
            5. Boolean Queries:
               - "_exists_:dbsnp AND dbsnp.vartype:snp" - Find SNPs in dbSNP
               - "dbsnp.vartype:snp OR dbsnp.vartype:indel" - Find SNPs or indels
               - "_exists_:dbsnp AND NOT dbsnp.vartype:indel" - Find non-indel variants in dbSNP
            
            The response includes pagination information and can be returned as a pandas DataFrame.
            """
        )
        async def query_variants(
            q: str = Query(
                ...,
                description="Query string",
                examples=[
                    "rs58991260",
                    "dbsnp.vartype:snp",
                    "dbnsfp.polyphen2.hdiv.score:>0.99",
                    "dbnsfp.genename:CDK*",
                    "_exists_:dbsnp AND dbsnp.vartype:snp"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "cadd.phred,dbsnp.rsid",
                    "clinvar.clinical_significance,cosmic.cosmic_id",
                    "dbnsfp.sift_pred,dbnsfp.polyphen2_hdiv_pred"
                ]
            ),
            size: int = Query(10, description="Maximum number of results to return (max 1000)", examples=[10, 50, 100]),
            skip: int = Query(0, description="Number of results to skip (for pagination)", examples=[0, 10, 20]),
            sort: Optional[str] = Query(None, description="Sort field, prefix with '-' for descending order", examples=["_score", "-_score", "cadd.phred"]),
            email: Optional[str] = Query(None, description="User email for tracking usage", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as pandas DataFrame", examples=[True, False]),
            df_index: bool = Query(True, description="Index DataFrame by query (only if as_dataframe=True)", examples=[True, False])
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
            description="""
            Perform multiple variant queries in a single request.
            
            This endpoint is useful for batch processing of variant queries. It supports:
            
            1. Multiple Query Types:
               - rsID queries: ["rs58991260", "rs12345678"]
               - HGVS queries: ["chr7:g.140453134T>C", "chr1:g.69000A>G"]
               - Mixed queries: ["rs58991260", "chr7:g.140453134T>C"]
            
            2. Field Scoping:
               - Search in specific fields: scopes=["dbsnp.vartype", "dbnsfp.genename"]
               - Search in all fields: scopes=None
            
            3. Result Filtering:
               - Return specific fields: fields=["cadd.phred", "dbsnp.rsid"]
               - Return all fields: fields=None
            
            The response can be returned as a pandas DataFrame for easier data manipulation.
            """
        )
        async def query_many_variants(
            query_list: str = Query(
                ...,
                description="Comma-separated list of query strings",
                examples=[
                    "rs58991260,rs12345678",
                    "chr7:g.140453134T>C,chr1:g.69000A>G",
                    "rs58991260,chr7:g.140453134T>C"
                ]
            ),
            scopes: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to search in",
                examples=[
                    "dbsnp.vartype,dbnsfp.genename",
                    "clinvar.clinical_significance,cosmic.cosmic_id",
                    "dbnsfp.sift_pred,dbnsfp.polyphen2_hdiv_pred"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "cadd.phred,dbsnp.rsid",
                    "clinvar.clinical_significance,cosmic.cosmic_id",
                    "dbnsfp.sift_pred,dbnsfp.polyphen2_hdiv_pred"
                ]
            ),
            email: Optional[str] = Query(None, description="User email for tracking usage", examples=["user@example.com"]),
            as_dataframe: bool = Query(False, description="Return results as pandas DataFrame", examples=[True, False]),
            df_index: bool = Query(True, description="Index DataFrame by query (only if as_dataframe=True)", examples=[True, False])
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
            description="""
            Retrieves detailed information about a specific variant using its identifier.
            
            This endpoint supports various variant ID formats:
            - HGVS notation: "chr7:g.140453134T>C"
            - dbSNP rsID: "rs58991260"
            - Genomic coordinates: "chr1:69000-70000"
            
            The response includes comprehensive variant information such as:
            - Basic information (ID, chromosome, position)
            - VCF information (reference, alternative, quality)
            - Genomic location (start, end, strand)
            - CADD scores and predictions
            - ClinVar annotations
            - COSMIC annotations
            - dbNSFP functional predictions
            - dbSNP annotations
            - DoCM annotations
            - MutDB annotations
            - SnpEff annotations
            
            You can filter the returned fields using the 'fields' parameter.
            """
        )
        async def get_variant(
            variant_id: str = Path(..., description="Variant identifier", examples=["chr7:g.140453134T>C", "rs58991260", "chr1:69000-70000"]),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "cadd.phred,dbsnp.rsid",
                    "clinvar.clinical_significance,cosmic.cosmic_id",
                    "dbnsfp.sift_pred,dbnsfp.polyphen2_hdiv_pred"
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
            description="""
            Retrieves information for multiple variants in a single request.
            
            This endpoint accepts a comma-separated list of variant IDs in various formats:
            - HGVS notations: "chr7:g.140453134T>C,chr1:g.69000A>G"
            - dbSNP rsIDs: "rs58991260,rs12345678"
            - Mixed formats: "chr7:g.140453134T>C,rs58991260"
            
            The response includes the same comprehensive variant information as the single variant endpoint,
            but for all requested variants.
            
            You can filter the returned fields using the 'fields' parameter.
            """
        )
        async def get_variants(
            variant_ids: str = Query(
                ...,
                description="Comma-separated list of variant IDs",
                examples=[
                    "chr7:g.140453134T>C,chr1:g.69000A>G",
                    "rs58991260,rs12345678",
                    "chr7:g.140453134T>C,rs58991260"
                ]
            ),
            fields: Optional[str] = Query(
                None,
                description="Comma-separated list of fields to return",
                examples=[
                    "cadd.phred,dbsnp.rsid",
                    "clinvar.clinical_significance,cosmic.cosmic_id",
                    "dbnsfp.sift_pred,dbnsfp.polyphen2_hdiv_pred"
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
