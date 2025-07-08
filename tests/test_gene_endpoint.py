import pytest
from biothings_mcp.server import BiothingsMCP
from biothings_mcp.biothings_api import GeneTools

@pytest.fixture
def mcp_server():
    """Fixture providing a BiothingsMCP server instance for testing."""
    return BiothingsMCP()

@pytest.fixture
def gene_tools(mcp_server):
    """Fixture providing GeneTools instance for testing."""
    return GeneTools(mcp_server, prefix="test_")

@pytest.mark.asyncio
async def test_get_gene(gene_tools):
    """Test the get_gene MCP tool.
    
    This test verifies that the tool correctly retrieves gene information
    for a given gene ID. It checks:
    1. The response contains the correct gene ID
    2. The response contains the correct gene symbol
    3. The response contains the correct gene name
    
    The test uses the CDK2 gene (ID: 1017) as an example.
    """
    result = await gene_tools.get_gene("1017")
    assert result.id == "1017"
    assert result.symbol == "CDK2"
    assert result.name == "cyclin dependent kinase 2"

@pytest.mark.asyncio
async def test_get_gene_with_fields(gene_tools):
    """Test the get_gene MCP tool with specific fields.
    
    This test verifies that the tool correctly filters the response to include
    only the requested fields. It checks:
    1. The response contains the requested fields (symbol and name)
    2. The response contains the required field (id)
    3. The field values are correct
    
    The test uses the CDK2 gene (ID: 1017) as an example and requests only
    the symbol and name fields.
    """
    result = await gene_tools.get_gene("1017", fields="symbol,name")
    # Check that we have the requested fields
    assert hasattr(result, "symbol")
    assert hasattr(result, "name")
    assert result.symbol == "CDK2"
    assert result.name == "cyclin dependent kinase 2"
    # Check that we have the required fields
    assert hasattr(result, "id")
    assert result.id == "1017"

@pytest.mark.asyncio
async def test_get_genes(gene_tools):
    """Test the get_genes MCP tool for multiple genes.
    
    This test verifies that the tool correctly retrieves information for
    multiple genes in a single request. It checks:
    1. The response contains the correct number of genes
    2. The response contains the correct gene IDs in the correct order
    
    The test uses two genes (CDK2 and CDK3) as examples.
    """
    result = await gene_tools.get_genes("1017,1018")
    assert len(result) == 2
    assert result[0].id == "1017"
    assert result[1].id == "1018"

@pytest.mark.asyncio
async def test_query_genes(gene_tools):
    """Test the query_genes MCP tool.
    
    This test verifies that the tool correctly performs a gene query and
    returns the expected results. It checks:
    1. The response contains a hits array
    2. The hits array contains at least one result
    3. The first hit contains the correct gene information
    
    The test queries for the CDK2 gene using its symbol.
    """
    result = await gene_tools.query_genes("symbol:CDK2", size=1)
    assert hasattr(result, "hits")
    assert len(result.hits) > 0
    hit = result.hits[0]
    assert hit.symbol == "CDK2"
    assert hit.name == "cyclin dependent kinase 2"

@pytest.mark.asyncio
async def test_query_many_genes(gene_tools):
    """Test the query_many_genes MCP tool.
    
    This test verifies that the tool correctly performs batch queries for
    multiple genes. It checks:
    1. The response contains the correct number of results
    2. The response contains results for all queried genes
    
    The test queries for two genes (CDK2 and BRCA1) using their symbols.
    """
    result = await gene_tools.query_many_genes("CDK2,BRCA1", scopes="symbol", size=1)
    assert len(result) == 2
    # Check that we got results for both queries
    symbols = {result_item.symbol for result_item in result}
    assert "CDK2" in symbols
    assert "BRCA1" in symbols

@pytest.mark.asyncio
async def test_gene_metadata(gene_tools):
    """Test the get_gene_metadata MCP tool.
    
    This test verifies that the tool correctly returns metadata about the
    gene database. It checks:
    1. The response is a MetadataResponse object
    2. The response contains a stats field
    3. The stats field is a dictionary
    4. The stats field contains a total field
    5. The total field is an integer
    
    The metadata includes information about the total number of genes in the
    database and other statistics.
    """
    result = await gene_tools.get_gene_metadata()
    assert hasattr(result, "stats")
    assert isinstance(result.stats, dict)
    assert "total" in result.stats
    assert isinstance(result.stats["total"], int)

@pytest.mark.asyncio
async def test_get_gene_ensembl_id(gene_tools):
    """Test the get_gene MCP tool with an Ensembl ID.

    This test verifies that the tool correctly retrieves gene information
    using an Ensembl gene ID. It checks:
    1. The response contains the correct gene ID (matches the queried Ensembl ID)
    2. Placeholder assertions for symbol and name (values need verification).

    The test uses the Ensembl ID ENSECAG00000002212 (likely from Cavia porcellus).
    """
    ensembl_id = "ENSECAG00000002212"
    result = await gene_tools.get_gene(ensembl_id)
    # Note: The actual ID field returned by MyGene.info for Ensembl might be different
    # (e.g., it might resolve to an Entrez ID or keep the Ensembl ID in a specific field).
    # We now check if the queried ensembl_id exists within the 'ensembl.gene' field.
    assert hasattr(result, "ensembl")
    assert "gene" in result.ensembl
    assert result.ensembl["gene"] == ensembl_id
    # TODO: Verify the expected symbol and name for this Ensembl ID and adjust assertions.
    # Check if standard fields are present, even if values aren't known
    assert hasattr(result, "symbol")
    assert hasattr(result, "name")
    assert hasattr(result, "taxid") 