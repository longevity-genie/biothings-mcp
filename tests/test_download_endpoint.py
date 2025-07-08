import pytest
import os
import tempfile
from pathlib import Path
from typing import Generator, Optional
from biothings_mcp.server import BiothingsMCP
from biothings_mcp.download_api import DownloadTools, get_entrez, run_pairwise_alignment, PairwiseAlignmentRequest, PairwiseAlignmentResponse, LocalFileResult

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
async def test_perform_pairwise_alignment_global(download_tools: DownloadTools) -> None:
    """Test the perform_pairwise_alignment MCP tool with global alignment.
    
    This test verifies that the tool correctly performs global pairwise
    sequence alignment and returns the expected results.
    """
    # Test data
    sequence1 = "GATTACA"
    sequence2 = "GCATGCU"
    
    request = PairwiseAlignmentRequest(
        sequence1=sequence1,
        sequence2=sequence2,
        match_score=2.0,
        mismatch_penalty=-1.0,
        open_gap_penalty=-1.0,
        extend_gap_penalty=-0.5,
        mode="global"
    )
    
    result: PairwiseAlignmentResponse = run_pairwise_alignment(request)
    
    # Check that we got a result with expected structure
    assert result is not None
    assert hasattr(result, "score")
    assert hasattr(result, "aligned_sequence1")
    assert hasattr(result, "aligned_sequence2")
    assert hasattr(result, "full_alignment_str")
    assert hasattr(result, "parameters_used")
    
    # Check that aligned sequences contain the original sequences
    assert len(result.aligned_sequence1) >= len(sequence1)
    assert len(result.aligned_sequence2) >= len(sequence2)
    
    # Check that parameters were stored correctly
    assert result.parameters_used["match_score"] == 2.0
    assert result.parameters_used["mismatch_penalty"] == -1.0
    assert result.parameters_used["mode"] == "global"

@pytest.mark.asyncio
async def test_perform_pairwise_alignment_local(download_tools: DownloadTools) -> None:
    """Test the perform_pairwise_alignment MCP tool with local alignment.
    
    This test verifies that the tool correctly performs local pairwise
    sequence alignment and returns the expected results.
    """
    # Test data - sequences with some similarity
    sequence1 = "ACGTACGTACGT"
    sequence2 = "CGTACGTA"
    
    request = PairwiseAlignmentRequest(
        sequence1=sequence1,
        sequence2=sequence2,
        match_score=1.0,
        mismatch_penalty=-1.0,
        open_gap_penalty=-0.5,
        extend_gap_penalty=-0.1,
        mode="local"
    )
    
    result: PairwiseAlignmentResponse = run_pairwise_alignment(request)
    
    # Check that we got a result with expected structure
    assert result is not None
    assert hasattr(result, "score")
    assert hasattr(result, "aligned_sequence1")
    assert hasattr(result, "aligned_sequence2")
    assert hasattr(result, "full_alignment_str")
    assert hasattr(result, "parameters_used")
    
    # Local alignment score should be positive for similar sequences
    assert result.score > 0
    
    # Check that parameters were stored correctly
    assert result.parameters_used["mode"] == "local"

@pytest.mark.asyncio
async def test_perform_pairwise_alignment_identical_sequences(download_tools: DownloadTools) -> None:
    """Test the perform_pairwise_alignment MCP tool with identical sequences.
    
    This test verifies that the tool correctly handles identical sequences
    and returns a perfect alignment score.
    """
    # Test data - identical sequences
    sequence = "ATCGATCG"
    
    request = PairwiseAlignmentRequest(
        sequence1=sequence,
        sequence2=sequence,
        match_score=1.0,
        mismatch_penalty=-1.0,
        mode="global"
    )
    
    result: PairwiseAlignmentResponse = run_pairwise_alignment(request)
    
    # Check that we got a result
    assert result is not None
    assert hasattr(result, "score")
    assert hasattr(result, "aligned_sequence1")
    assert hasattr(result, "aligned_sequence2")
    
    # Identical sequences should have perfect alignment
    assert result.aligned_sequence1 == sequence
    assert result.aligned_sequence2 == sequence
    # Score should be positive (number of matches * match_score)
    assert result.score > 0

@pytest.mark.asyncio
async def test_perform_pairwise_alignment_local_save_to_file(download_tools: DownloadTools) -> None:
    """Test the perform_pairwise_alignment_local functionality.
    
    This test verifies that the tool correctly performs pairwise alignment
    and saves the results to a local file.
    """
    # Test data
    sequence1 = "GATTACA"
    sequence2 = "GCATGCU"
    
    request = PairwiseAlignmentRequest(
        sequence1=sequence1,
        sequence2=sequence2,
        match_score=2.0,
        mismatch_penalty=-1.0,
        open_gap_penalty=-1.0,
        extend_gap_penalty=-0.5,
        mode="global"
    )
    
    # Get the alignment result
    alignment_result: PairwiseAlignmentResponse = run_pairwise_alignment(request)
    
    # Create detailed alignment output (similar to what the local tool would create)
    alignment_content = f"""Pairwise Alignment Results
=============================

Parameters:
- Match Score: {alignment_result.parameters_used['match_score']}
- Mismatch Penalty: {alignment_result.parameters_used['mismatch_penalty']}
- Open Gap Penalty: {alignment_result.parameters_used['open_gap_penalty']}
- Extend Gap Penalty: {alignment_result.parameters_used['extend_gap_penalty']}
- Mode: {alignment_result.parameters_used['mode']}

Alignment Score: {alignment_result.score}

Aligned Sequences:
{alignment_result.full_alignment_str}

Sequence 1 (aligned): {alignment_result.aligned_sequence1}
Sequence 2 (aligned): {alignment_result.aligned_sequence2}
"""
    
    # Test saving to local file
    result: LocalFileResult = download_tools._save_to_local_file(
        data=alignment_content,
        format_type="alignment",
        output_path=None,
        default_prefix="pairwise_alignment"
    )
    
    # Check that the file was created successfully
    assert result["success"] is True
    assert result["format"] == "alignment"
    assert result["path"] is not None
    
    # Check that the file exists and contains the expected content
    file_path = Path(result["path"])
    assert file_path.exists()
    assert file_path.suffix == ".aln"
    
    with open(file_path, 'r') as f:
        saved_content = f.read()
    
    assert saved_content == alignment_content
    assert "Pairwise Alignment Results" in saved_content
    assert f"Alignment Score: {alignment_result.score}" in saved_content
    assert "GATTACA" in saved_content
    assert "GCATGCU" in saved_content

@pytest.mark.asyncio
async def test_output_directory_creation(download_tools: DownloadTools) -> None:
    """Test that output directory is created correctly."""
    # Check that output directory exists
    assert download_tools.output_dir.exists()
    assert download_tools.output_dir.is_dir()

@pytest.mark.asyncio
async def test_save_to_local_file_json_format(download_tools: DownloadTools) -> None:
    """Test saving data in JSON format."""
    test_data = {
        "gene_id": "ENSG00000141510",
        "gene_name": "TP53",
        "description": "tumor protein p53"
    }
    
    result: LocalFileResult = download_tools._save_to_local_file(
        data=test_data,
        format_type="json",
        output_path="test_gene_data",
        default_prefix="test"
    )
    
    # Check that the file was created successfully
    assert result["success"] is True
    assert result["format"] == "json"
    assert result["path"] is not None
    
    # Check that the file exists and contains the expected content
    file_path = Path(result["path"])
    assert file_path.exists()
    assert file_path.suffix == ".json"
    
    import json
    with open(file_path, 'r') as f:
        saved_data = json.load(f)
    
    assert saved_data == test_data

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
