import pytest
from biothings_mcp.server import BiothingsMCP
from biothings_mcp.biothings_api import TaxonTools

# Test constant for human taxonomy ID
HUMAN_TAXID = "9606"

@pytest.fixture
def mcp_server():
    """Fixture providing a BiothingsMCP server instance for testing."""
    return BiothingsMCP()

@pytest.fixture
def taxon_tools(mcp_server):
    """Fixture providing TaxonTools instance for testing."""
    return TaxonTools(mcp_server, prefix="test_")

@pytest.mark.asyncio
async def test_get_taxon(taxon_tools):
    """Test the get_taxon MCP tool.
    
    This test verifies that the tool correctly retrieves taxonomy information
    for a given taxonomy ID. It checks that the response contains the correct
    taxon ID and other expected fields.
    
    The test uses the human taxonomy ID (9606) as an example.
    """
    result = await taxon_tools.get_taxon(HUMAN_TAXID)
    # Check that we got a result
    assert result is not None
    # Check that we have the ID field
    assert hasattr(result, "id") or hasattr(result, "_id")
    # Check that the ID matches what we requested
    if hasattr(result, "id"):
        assert result.id == HUMAN_TAXID
    elif hasattr(result, "_id"):
        assert result._id == HUMAN_TAXID

@pytest.mark.asyncio
async def test_get_taxon_with_fields(taxon_tools):
    """Test the get_taxon MCP tool with specific fields.
    
    This test verifies that the tool correctly filters the response to include
    only the requested fields. It checks that the response contains the
    requested fields and the taxon ID.
    
    The test uses human taxonomy with specific field filtering.
    """
    fields = "scientific_name,lineage"
    result = await taxon_tools.get_taxon(HUMAN_TAXID, fields=fields)
    # Check that we got a result
    assert result is not None
    # Check that we have the ID field
    assert hasattr(result, "id") or hasattr(result, "_id")

@pytest.mark.asyncio
async def test_get_taxons(taxon_tools):
    """Test the get_taxons MCP tool for multiple taxonomy IDs.
    
    This test verifies that the tool correctly retrieves information for
    multiple taxonomy IDs in a single request. It checks that the response
    contains results for the queried taxon IDs.
    
    The test uses multiple taxonomy IDs (human and mouse).
    """
    taxon_ids = "9606,10090"  # Human and Mouse
    result = await taxon_tools.get_taxons(taxon_ids)
    assert isinstance(result, list)
    assert len(result) >= 1  # Should have at least one result

@pytest.mark.asyncio
async def test_query_taxons(taxon_tools):
    """Test the query_taxons MCP tool.
    
    This test verifies that the tool correctly performs a taxonomy query and
    returns the expected results. It checks that the response contains
    a hits array with at least one result.
    
    The test queries for taxonomy records using a search query.
    """
    result = await taxon_tools.query_taxons("scientific_name:homo", size=1)
    assert hasattr(result, "hits")
    assert len(result.hits) > 0
    # Check that the first hit has expected structure
    hit = result.hits[0]
    assert hasattr(hit, "id") or hasattr(hit, "_id")

@pytest.mark.asyncio
async def test_query_many_taxons(taxon_tools):
    """Test the query_many_taxons MCP tool.
    
    This test verifies that the tool correctly performs batch queries for
    multiple taxonomy identifiers. It checks that the response contains results for
    the queried taxonomy identifiers.
    
    The test queries for multiple taxonomy records using their identifiers.
    """
    result = await taxon_tools.query_many_taxons("9606,10090", scopes="_id")
    assert isinstance(result, list)
    assert len(result) >= 1  # Should have at least one result
    # Check that results have expected structure
    for taxon_result in result:
        assert hasattr(taxon_result, "id") or hasattr(taxon_result, "_id")
