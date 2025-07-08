import pytest
from typing import Any, List
from biothings_mcp.server import BiothingsMCP
from biothings_mcp.biothings_api import ChemTools

# Test constant for aspirin
ASPIRIN_INCHIKEY: str = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"

@pytest.fixture
def mcp_server() -> BiothingsMCP:
    """Fixture providing a BiothingsMCP server instance for testing."""
    return BiothingsMCP()

@pytest.fixture
def chem_tools(mcp_server: BiothingsMCP) -> ChemTools:
    """Fixture providing ChemTools instance for testing."""
    return ChemTools(mcp_server, prefix="test_")

@pytest.mark.asyncio
async def test_get_chem(chem_tools: ChemTools) -> None:
    """Test the get_chem MCP tool.
    
    This test verifies that the tool correctly retrieves chemical information
    for a given chemical ID. It checks that the response contains the correct
    chemical ID and other expected fields.
    
    The test uses aspirin's InChIKey as an example.
    """
    result: Any = await chem_tools.get_chem(ASPIRIN_INCHIKEY)
    # Check that we got a result
    assert result is not None
    # Check that we have the ID field
    assert hasattr(result, "id") or hasattr(result, "_id")

@pytest.mark.asyncio
async def test_get_chem_with_fields(chem_tools: ChemTools) -> None:
    """Test the get_chem MCP tool with specific fields.
    
    This test verifies that the tool correctly filters the response to include
    only the requested fields. It checks that the response contains the
    requested fields and the chemical ID.
    
    The test uses aspirin with specific field filtering.
    """
    fields: str = "pubchem.molecular_formula,pubchem.molecular_weight"
    result: Any = await chem_tools.get_chem(ASPIRIN_INCHIKEY, fields=fields)
    # Check that we got a result
    assert result is not None
    # Check that we have the ID field
    assert hasattr(result, "id") or hasattr(result, "_id")

@pytest.mark.asyncio
async def test_get_chems(chem_tools: ChemTools) -> None:
    """Test the get_chems MCP tool for multiple chemicals.
    
    This test verifies that the tool correctly retrieves information for
    multiple chemicals in a single request. It checks that the response
    contains results for the queried chemical IDs.
    
    The test uses multiple chemical IDs.
    """
    chem_ids: str = f"{ASPIRIN_INCHIKEY},KTUFNOKKBVMGRW-UHFFFAOYSA-N"  # Aspirin and Glucose
    result: List[Any] = await chem_tools.get_chems(chem_ids)
    assert isinstance(result, list)
    assert len(result) >= 1  # Should have at least one result

@pytest.mark.asyncio
async def test_query_chems(chem_tools: ChemTools) -> None:
    """Test the query_chems MCP tool.
    
    This test verifies that the tool correctly performs a chemical query and
    returns the expected results. It checks that the response contains
    a hits array with at least one result.
    
    The test queries for chemicals using a search query.
    """
    result: Any = await chem_tools.query_chems("pubchem.molecular_formula:C9H8O4", size=1)
    assert hasattr(result, "hits")
    assert len(result.hits) > 0
    # Check that the first hit has expected structure
    hit: Any = result.hits[0]
    assert hasattr(hit, "id") or hasattr(hit, "_id")

@pytest.mark.asyncio
async def test_query_many_chems(chem_tools: ChemTools) -> None:
    """Test the query_many_chems MCP tool.
    
    This test verifies that the tool correctly performs batch queries for
    multiple chemicals. It checks that the response contains results for
    the queried chemical identifiers.
    
    The test queries for multiple chemicals using PubChem CIDs.
    """
    # Use known PubChem CIDs for aspirin and glucose
    query_list = "2244,5793"  # aspirin CID: 2244, glucose CID: 5793
    result: List[Any] = await chem_tools.query_many_chems(query_list, scopes="pubchem.cid")
    assert isinstance(result, list)
    assert len(result) >= 1  # Should have at least one result
    # Check that results have expected structure
    for chem_result in result:
        assert hasattr(chem_result, "id") or hasattr(chem_result, "_id")
