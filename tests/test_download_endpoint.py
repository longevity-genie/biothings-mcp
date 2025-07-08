import pytest
import os
import tempfile
from pathlib import Path
from typing import Generator, Optional
from biothings_mcp.server import BiothingsMCP
from biothings_mcp.download_api import DownloadTools, get_entrez, LocalFileResult

@pytest.fixture
def temp_output_dir() -> Generator[str, None, None]:
    """Fixture providing a temporary output directory for testing."""
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir

@pytest.fixture
def mcp_server(temp_output_dir: str) -> BiothingsMCP:
    """Fixture providing a BiothingsMCP server instance for testing."""
    return BiothingsMCP(output_dir=temp_output_dir)

@pytest.fixture
def download_tools(mcp_server: BiothingsMCP) -> DownloadTools:
    """Fixture providing DownloadTools instance for testing."""
    return DownloadTools(mcp_server, prefix="test_", output_dir=mcp_server.output_dir)

@pytest.mark.asyncio
async def test_download_entrez_data_fasta(download_tools: DownloadTools) -> None:
    """Test the download_entrez_data MCP tool with FASTA format.
    
    This test verifies that the tool correctly downloads sequence data
    from NCBI Entrez in FASTA format. It checks that the response
    contains valid FASTA content.
    
    Note: This test requires ENTREZ_EMAIL environment variable to be set.
    """
    # Skip if no ENTREZ_EMAIL is set
    if not os.getenv("ENTREZ_EMAIL"):
        pytest.skip("ENTREZ_EMAIL environment variable not set")
    
    # Test data
    ids = ["NM_000546.6"]  # Human TP53 mRNA
    db = "nucleotide"
    reftype = "fasta"
    
    result: str = get_entrez(ids=ids, db=db, reftype=reftype)
    
    # Check that we got a result
    assert result is not None
    assert isinstance(result, str)
    # Basic FASTA format check
    assert result.startswith(">")
    assert "NM_000546.6" in result

@pytest.mark.asyncio
async def test_download_entrez_data_genbank(download_tools: DownloadTools) -> None:
    """Test the download_entrez_data MCP tool with GenBank format.
    
    This test verifies that the tool correctly downloads sequence data
    from NCBI Entrez in GenBank format. It checks that the response
    contains valid GenBank content.
    
    Note: This test requires ENTREZ_EMAIL environment variable to be set.
    """
    # Skip if no ENTREZ_EMAIL is set
    if not os.getenv("ENTREZ_EMAIL"):
        pytest.skip("ENTREZ_EMAIL environment variable not set")
    
    # Test data
    ids = ["NM_000546.6"]  # Human TP53 mRNA
    db = "nucleotide"
    reftype = "gb"
    
    result: str = get_entrez(ids=ids, db=db, reftype=reftype)
    
    # Check that we got a result
    assert result is not None
    assert isinstance(result, str)
    # Basic GenBank format check
    assert "LOCUS" in result
    assert "NM_000546.6" in result

@pytest.mark.asyncio
async def test_download_entrez_multiple_ids(download_tools: DownloadTools) -> None:
    """Test the download_entrez_data MCP tool with multiple IDs.
    
    This test verifies that the tool correctly downloads multiple sequence records
    in a single request.
    
    Note: This test requires ENTREZ_EMAIL environment variable to be set.
    """
    # Skip if no ENTREZ_EMAIL is set
    if not os.getenv("ENTREZ_EMAIL"):
        pytest.skip("ENTREZ_EMAIL environment variable not set")
    
    # Test data - multiple sequences
    ids = ["NM_000546.6", "NM_001126112.3"]  # Two different sequences
    db = "nucleotide"
    reftype = "fasta"
    
    result: str = get_entrez(ids=ids, db=db, reftype=reftype)
    
    # Check that we got a result
    assert result is not None
    assert isinstance(result, str)
    # Should contain both sequences
    assert "NM_000546.6" in result
    assert "NM_001126112.3" in result
    # Should have multiple FASTA headers
    assert result.count(">") >= 2

@pytest.mark.asyncio
async def test_download_entrez_data_local_fasta(download_tools: DownloadTools) -> None:
    """Test the download_entrez_data_local MCP tool with FASTA format.
    
    This test verifies that the tool correctly downloads sequence data
    from NCBI Entrez and saves it to a local file.
    
    Note: This test requires ENTREZ_EMAIL environment variable to be set.
    """
    # Skip if no ENTREZ_EMAIL is set
    if not os.getenv("ENTREZ_EMAIL"):
        pytest.skip("ENTREZ_EMAIL environment variable not set")
    
    # Test data
    ids = ["NM_000546.6"]  # Human TP53 mRNA
    db = "nucleotide"
    reftype = "fasta"
    
    # Get the data first
    content: str = get_entrez(ids=ids, db=db, reftype=reftype)
    
    # Test the save to local file functionality
    result: LocalFileResult = download_tools._save_to_local_file(
        data=content,
        format_type=reftype,
        output_path=None,
        default_prefix=f"entrez_{db}"
    )
    
    # Check that the file was created successfully
    assert result["success"] is True
    assert result["format"] == reftype
    assert result["path"] is not None
    
    # Check that the file exists and contains the expected content
    file_path = Path(result["path"])
    assert file_path.exists()
    assert file_path.suffix == ".fasta"
    
    with open(file_path, 'r') as f:
        saved_content = f.read()
    
    assert saved_content == content
    assert saved_content.startswith(">")
    assert "NM_000546.6" in saved_content

@pytest.mark.asyncio
async def test_download_entrez_data_local_genbank(download_tools: DownloadTools) -> None:
    """Test the download_entrez_data_local MCP tool with GenBank format.
    
    This test verifies that the tool correctly downloads sequence data
    from NCBI Entrez in GenBank format and saves it to a local file.
    
    Note: This test requires ENTREZ_EMAIL environment variable to be set.
    """
    # Skip if no ENTREZ_EMAIL is set
    if not os.getenv("ENTREZ_EMAIL"):
        pytest.skip("ENTREZ_EMAIL environment variable not set")
    
    # Test data
    ids = ["NM_000546.6"]  # Human TP53 mRNA
    db = "nucleotide"
    reftype = "gb"
    
    # Get the data first
    content: str = get_entrez(ids=ids, db=db, reftype=reftype)
    
    # Test the save to local file functionality
    result: LocalFileResult = download_tools._save_to_local_file(
        data=content,
        format_type=reftype,
        output_path=None,
        default_prefix=f"entrez_{db}"
    )
    
    # Check that the file was created successfully
    assert result["success"] is True
    assert result["format"] == reftype
    assert result["path"] is not None
    
    # Check that the file exists and contains the expected content
    file_path = Path(result["path"])
    assert file_path.exists()
    assert file_path.suffix == ".gb"
    
    with open(file_path, 'r') as f:
        saved_content = f.read()
    
    assert saved_content == content
    assert "LOCUS" in saved_content
    assert "NM_000546.6" in saved_content

@pytest.mark.asyncio
async def test_download_entrez_data_local_custom_path(download_tools: DownloadTools) -> None:
    """Test the download_entrez_data_local MCP tool with custom output path.
    
    This test verifies that the tool correctly saves data to a custom output path.
    
    Note: This test requires ENTREZ_EMAIL environment variable to be set.
    """
    # Skip if no ENTREZ_EMAIL is set
    if not os.getenv("ENTREZ_EMAIL"):
        pytest.skip("ENTREZ_EMAIL environment variable not set")
    
    # Test data
    ids = ["NM_000546.6"]  # Human TP53 mRNA
    db = "nucleotide"
    reftype = "fasta"
    custom_filename = "custom_tp53_sequence"
    
    # Get the data first
    content: str = get_entrez(ids=ids, db=db, reftype=reftype)
    
    # Test the save to local file functionality with custom path
    result: LocalFileResult = download_tools._save_to_local_file(
        data=content,
        format_type=reftype,
        output_path=custom_filename,
        default_prefix=f"entrez_{db}"
    )
    
    # Check that the file was created successfully
    assert result["success"] is True
    assert result["format"] == reftype
    assert result["path"] is not None
    
    # Check that the file exists with the custom name
    file_path = Path(result["path"])
    assert file_path.exists()
    assert file_path.name == f"{custom_filename}.fasta"
    
    with open(file_path, 'r') as f:
        saved_content = f.read()
    
    assert saved_content == content


@pytest.mark.asyncio
async def test_save_to_local_file_error_handling(download_tools: DownloadTools) -> None:
    """Test error handling in _save_to_local_file method."""
    # Test with invalid path (trying to write to a directory that doesn't exist and can't be created)
    invalid_path = "/root/nonexistent/path/test_file"
    
    result: LocalFileResult = download_tools._save_to_local_file(
        data="test content",
        format_type="txt",
        output_path=invalid_path,
        default_prefix="test"
    )
    
    # Should handle the error gracefully
    assert result["success"] is False
    assert result["path"] is None
    assert "error" in result
    assert result["format"] == "txt"
