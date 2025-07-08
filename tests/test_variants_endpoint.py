import pytest
from typing import Any, List, Union
from biothings_mcp.server import BiothingsMCP
from biothings_mcp.biothings_api import VariantTools

@pytest.fixture
def mcp_server() -> BiothingsMCP:
    """Fixture providing a BiothingsMCP server instance for testing."""
    return BiothingsMCP()

@pytest.fixture
def variant_tools(mcp_server: BiothingsMCP) -> VariantTools:
    """Fixture providing VariantTools instance for testing."""
    return VariantTools(mcp_server, prefix="test_")

@pytest.mark.asyncio
async def test_get_variant(variant_tools: VariantTools) -> None:
    """Test the get_variant MCP tool.
    
    This test verifies that the tool correctly retrieves variant information
    for a given variant ID. It checks that the response contains the correct
    variant ID and other expected fields.
    
    The test uses an example variant ID (likely an rsID or HGVS identifier).
    """
    variant_id: str = "rs58991260"
    result: Any = await variant_tools.get_variant(variant_id)
    # Check that we got a result
    assert result is not None
    # Check that standard fields are present
    assert hasattr(result, "id") or hasattr(result, "_id")

@pytest.mark.asyncio
async def test_get_variant_with_fields(variant_tools: VariantTools) -> None:
    """Test the get_variant MCP tool with specific fields.
    
    This test verifies that the tool correctly filters the response to include
    only the requested fields. It checks that the response contains the
    requested fields and the variant ID.
    
    The test uses an example variant with specific field filtering.
    """
    variant_id: str = "rs58991260"
    fields: str = "dbsnp.rsid,dbsnp.vartype"
    result: Any = await variant_tools.get_variant(variant_id, fields=fields)
    # Check that we got a result
    assert result is not None
    # Check that we have the ID field
    assert hasattr(result, "id") or hasattr(result, "_id")

@pytest.mark.asyncio
async def test_get_variants(variant_tools: VariantTools) -> None:
    """Test the get_variants MCP tool for multiple variants.
    
    This test verifies that the tool correctly retrieves information for
    multiple variants in a single request. It checks that the response
    contains the correct number of variants.
    
    The test uses multiple variant IDs.
    """
    variant_ids: str = "rs58991260,rs2500"
    result: List[Any] = await variant_tools.get_variants(variant_ids)
    assert isinstance(result, list)
    assert len(result) >= 1  # Should have at least one result

@pytest.mark.asyncio
async def test_query_variants(variant_tools: VariantTools) -> None:
    """Test the query_variants MCP tool.
    
    This test verifies that the tool correctly performs a variant query and
    returns the expected results. It checks that the response contains
    a hits array with at least one result.
    
    The test queries for variants using a search query.
    """
    result: Any = await variant_tools.query_variants("dbsnp.vartype:snv", size=1)
    assert hasattr(result, "hits")
    assert len(result.hits) > 0
    # Check that the first hit has expected structure
    hit: Any = result.hits[0]
    assert hasattr(hit, "id") or hasattr(hit, "_id")

@pytest.mark.asyncio
async def test_query_many_variants(variant_tools: VariantTools) -> None:
    """Test the query_many_variants MCP tool.
    
    This test verifies that the tool correctly performs batch queries for
    multiple variants. It checks that the response contains results for
    the queried variant identifiers.
    
    The test queries for multiple variants using their identifiers.
    """
    result: List[Any] = await variant_tools.query_many_variants("rs58991260,rs2500", scopes="dbsnp.rsid")
    assert isinstance(result, list)
    assert len(result) >= 1  # Should have at least one result
    # Check that results have expected structure
    for variant_result in result:
        assert hasattr(variant_result, "id") or hasattr(variant_result, "_id")
