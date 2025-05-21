import pytest
from fastapi.testclient import TestClient
from biothings_mcp.server import create_app
from pathlib import Path
import logging
import os

# Configure logger for this test module
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG) # Or logging.INFO for less verbosity

# It's good practice to ensure ENTREZ_EMAIL is set for tests that might call NCBI
# For testing purposes, a placeholder can be used if not set, but NCBI might block it.
# Consider mocking Entrez calls for robust CI/CD if actual calls are problematic.
if "ENTREZ_EMAIL" not in os.environ:
    os.environ["ENTREZ_EMAIL"] = "test_user@example.com"
    logger.warning(f"ENTREZ_EMAIL not set, using placeholder: {os.environ['ENTREZ_EMAIL']}")


@pytest.fixture
def client():
    """Fixture providing a FastAPI test client."""
    project_root = Path(__file__).resolve().parents[1]  # Project root is one level up from tests dir
    log_dir = project_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)  # Ensure the directory exists

    app = create_app()
    return TestClient(app)


def test_entrez_download_fasta_nucleotide(client):
    """Test downloading a nucleotide sequence in FASTA format from Entrez."""
    # Using a known human mRNA sequence for TP53
    payload = {
        "ids": ["NM_000546.6"],
        "db": "nucleotide",
        "reftype": "fasta"
    }
    response = client.post("/download/entrez", json=payload)
    assert response.status_code == 200
    content = response.json()
    assert content.startswith(">NM_000546.6 Homo sapiens tumor protein p53 (TP53), transcript variant 1, mRNA")
    # Basic check for sequence content (can be more specific if needed)
    assert "GATTACA" in content or "ACGT" in content # Generic DNA subsequence check
    logger.info(f"Entrez FASTA download successful for NM_000546.6. Length: {len(content)}")

def test_entrez_download_genbank_protein(client):
    """Test downloading a protein sequence in GenBank format from Entrez."""
    # Using a known human p53 protein sequence
    payload = {
        "ids": ["NP_000537.3"],
        "db": "protein",
        "reftype": "gb"  # GenBank format
    }
    response = client.post("/download/entrez", json=payload)
    assert response.status_code == 200
    content = response.json()
    assert "LOCUS       NP_000537" in content
    assert "DEFINITION  cellular tumor antigen p53 isoform a [Homo sapiens]." in content
    assert "ORGANISM  Homo sapiens" in content
    logger.info(f"Entrez GenBank download successful for NP_000537.3. Length: {len(content)}")

    # Verify the actual sequence
    expected_sequence = (
        "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAA"
        "PPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKT"
        "CPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRN"
        "TFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGR"
        "DRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALEL"
        "KDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD"
    )
    # Helper to extract sequence from GenBank format
    def extract_sequence_from_genbank(gb_content):
        sequence_lines = []
        in_sequence_section = False
        for line in gb_content.splitlines():
            if line.startswith("ORIGIN"):
                in_sequence_section = True
                continue
            if in_sequence_section:
                if line.startswith("//"): # End of record
                    break
                # Remove line numbers and spaces, then concatenate sequence parts
                parts = line.split()
                if len(parts) > 1: # Ensure there's sequence data beyond the number
                    sequence_lines.append("".join(parts[1:]))
        return "".join(sequence_lines).upper() # NCBI sequences are typically uppercase

    actual_sequence = extract_sequence_from_genbank(content)
    assert actual_sequence == expected_sequence, "The downloaded protein sequence does not match the expected sequence."
    logger.info(f"Entrez GenBank sequence verification successful for NP_000537.3.")

def test_entrez_download_invalid_id(client):
    """Test Entrez download with an invalid ID."""
    payload = {
        "ids": ["INVALID_ID_123"],
        "db": "nucleotide",
        "reftype": "fasta"
    }
    response = client.post("/download/entrez", json=payload)
    # NCBI often returns 400 for bad requests like invalid ID format
    # or if ID does not exist, the efetch will return an empty result which might not be an error from client perspective
    # but here the API should probably propagate an error if Entrez.efetch itself fails or returns empty unexpectedly.
    # Based on current implementation, Entrez.efetch might return empty string for invalid ID without raising HTTPError directly.
    # The API wrapper in `download_api.py` might need adjustment to explicitly check for empty results if that's an error condition.
    # For now, let's assume NCBI might return 200 with empty or minimal content for non-existent specific IDs if no other error occurs.
    # If the ID format is grossly invalid, a 400 is more likely from NCBI directly.
    # Let's test for 200 and empty content, or a 400/500 if the API handles it as an error.
    if response.status_code == 200:
        assert response.json() == "", "Expected empty response for invalid ID if status is 200"
        logger.warning("Entrez download with invalid ID returned 200 with empty content.")
    else:
        assert response.status_code in [400, 404, 500] # More robust check for various error codes
        logger.info(f"Entrez download with invalid ID correctly failed with status: {response.status_code}")

def test_entrez_download_invalid_db(client):
    """Test Entrez download with an invalid database name."""
    payload = {
        "ids": ["NM_000546.6"],
        "db": "invalid_database_name",
        "reftype": "fasta"
    }
    response = client.post("/download/entrez", json=payload)
    # Expecting a 400 or 500 due to server-side validation or error from Entrez
    assert response.status_code >= 400
    data = response.json()
    assert "detail" in data
    logger.info(f"Entrez download with invalid DB correctly failed with status: {response.status_code}, details: {data.get('detail')}")

def test_pairwise_alignment_global_valid(client):
    """Test global pairwise alignment with valid sequences."""
    payload = {
        "sequence1": "GATTACA",
        "sequence2": "GCATGCU", # Note: U is not standard DNA, Biopython handles it
        "match_score": 2.0,
        "mismatch_penalty": -1.0,
        "open_gap_penalty": -1.0,
        "extend_gap_penalty": -0.5,
        "mode": "global"
    }
    response = client.post("/tools/align/pairwise", json=payload)
    assert response.status_code == 200, f"Expected 200 OK but got {response.status_code}. Response: {response.text}"
    data = response.json()

    assert "score" in data
    assert isinstance(data["score"], float)
    assert "aligned_sequence1" in data
    assert isinstance(data["aligned_sequence1"], str)
    assert "aligned_sequence2" in data
    assert isinstance(data["aligned_sequence2"], str)
    assert "full_alignment_str" in data
    assert isinstance(data["full_alignment_str"], str)
    assert "parameters_used" in data
    assert isinstance(data["parameters_used"], dict)
    assert data["parameters_used"]["sequence1"] == payload["sequence1"] # Check if params are echoed back

    # Basic check that aligned sequences are not empty and have similar lengths (can be refined)
    assert len(data["aligned_sequence1"]) > 0
    assert len(data["aligned_sequence2"]) > 0
    assert len(data["aligned_sequence1"]) == len(data["aligned_sequence2"])

    # Log the actual alignment for review
    logger.info(f"Global alignment successful. Score: {data['score']}")
    logger.info(f"Alignment:\n{data['full_alignment_str']}")
    # Example of a more specific check (this would require knowing the exact expected alignment)
    # expected_aligned1 = "GATTACA" 
    # expected_aligned2 = "GCA-TGC" # This is an example, actual alignment depends on scoring
    # assert data["aligned_sequence1"] == expected_aligned1
    # assert data["aligned_sequence2"] == expected_aligned2
    # assert data["score"] ==  -0.5 # Example score, needs to be calculated or known

def test_pairwise_alignment_local_valid(client):
    """Test local pairwise alignment with valid sequences."""
    payload = {
        "sequence1": "AGCTAGCTAGCT",
        "sequence2": "GCTAGC",
        "match_score": 5.0,
        "mismatch_penalty": -4.0,
        "open_gap_penalty": -10.0,
        "extend_gap_penalty": -0.5,
        "mode": "local"
    }
    response = client.post("/tools/align/pairwise", json=payload)
    assert response.status_code == 200
    data = response.json()
    assert "score" in data
    assert isinstance(data["score"], float)
    assert "aligned_sequence1" in data
    assert isinstance(data["aligned_sequence1"], str)
    assert "aligned_sequence2" in data
    assert isinstance(data["aligned_sequence2"], str)
    assert "full_alignment_str" in data
    assert isinstance(data["full_alignment_str"], str)
    assert "parameters_used" in data
    assert isinstance(data["parameters_used"], dict)
    logger.info(f"Local alignment request resulted in 200 as expected.")

def test_pairwise_alignment_empty_sequence(client):
    """Test pairwise alignment with one empty sequence."""
    payload = {
        "sequence1": "",
        "sequence2": "GCATGCU",
        "mode": "global"
    }
    response = client.post("/tools/align/pairwise", json=payload)
    # Biopython might produce a result (e.g. aligning to all gaps) or the API might raise an error.
    # The current implementation of `run_pairwise_alignment` would likely lead to a ValueError from Biopython if sequences are empty or invalid.
    # This ValueError is caught and results in a 400 error.
    assert response.status_code == 400
    # Check for an appropriate error message if the API provides one
    data = response.json()
    assert "detail" in data
    assert isinstance(data["detail"], str) # Expect a string detail
    assert "sequence has zero length" in data["detail"] # Check for specific error message
    logger.info(f"Pairwise alignment with empty sequence correctly failed with status 400 and detail: {data['detail']}")

def test_pairwise_alignment_invalid_params(client):
    """Test pairwise alignment with invalid scoring parameters (e.g., positive mismatch penalty if not intended)."""
    payload = {
        "sequence1": "GATTACA",
        "sequence2": "GCATGCU",
        "match_score": 1.0,
        "mismatch_penalty": 1.0,  # Positive mismatch penalty, usually not desired/supported explicitly by some algos
        "open_gap_penalty": 0.5, # Positive open gap penalty
        "extend_gap_penalty": 0.1, # Positive extend gap penalty
        "mode": "global"
    }
    # The Pydantic model has constraints like "Should be negative or zero."
    # However, those are descriptions, not hard Pydantic validation constraints like `lt=0`.
    # Biopython itself might accept these values but produce unusual results, or it might raise an error for some combinations.
    # The API itself does not add further validation on these values beyond Pydantic's type checks.
    # Let's assume Biopython handles it, and we check for a 200, but the results might be "unphysical".
    response = client.post("/tools/align/pairwise", json=payload)
    # If Biopython raises an error due to parameters, it would be a 400 or 500 from the API.
    # For this test, we'll assume it processes and returns a result.
    # A more robust test might check specific Biopython behavior for such params if known.
    assert response.status_code == 200 # or 400 if Biopython/Pydantic eventually rejects it
