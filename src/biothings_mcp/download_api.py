import os
from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from pydantic import BaseModel, Field, ConfigDict
from pathlib import Path
from typing import Literal, List, Dict, Optional
from fastapi import HTTPException
from urllib.error import HTTPError


DB_LITERAL = Literal[
    "pubmed", "protein", "nuccore", "ipg", "nucleotide", "structure", "genome",
    "annotinfo", "assembly", "bioproject", "biosample", "blastdbinfo", "books",
    "cdd", "clinvar", "gap", "gapplus", "grasp", "dbvar", "gene", "gds",
    "geoprofiles", "medgen", "mesh", "nlmcatalog", "omim", "orgtrack", "pmc",
    "proteinclusters", "pcassay", "protfam", "pccompound", "pcsubstance",
    "seqannot", "snp", "sra", "taxonomy", "biocollections", "gtr"
]


class EntrezDownloadRequest(BaseModel):
    ids: List[str]
    db: DB_LITERAL
    reftype: Literal["fasta", "gb"]

    model_config = {
        "json_schema_extra": {
            "example": {
                "ids": ["NM_000546.6"],
                "db": "nucleotide",
                "reftype": "fasta",
            }
        }
    }


class PairwiseAlignmentRequest(BaseModel):
    sequence1: str = Field(..., description="First sequence for alignment.")
    sequence2: str = Field(..., description="Second sequence for alignment.")
    match_score: float = Field(1.0, description="Score for a match.")
    mismatch_penalty: float = Field(-1.0, description="Penalty for a mismatch. Should be negative or zero.")
    open_gap_penalty: float = Field(-0.5, description="Penalty for opening a gap. Should be negative or zero.")
    extend_gap_penalty: float = Field(-0.1, description="Penalty for extending a gap. Should be negative or zero.")
    mode: Literal["global", "local"] = Field("global", description="Alignment mode: 'global' or 'local'.")

    model_config = ConfigDict(
        json_schema_extra = {
            "example": {
                "sequence1": "GATTACA",
                "sequence2": "GCATGCU",
                "match_score": 2.0,
                "mismatch_penalty": -1.0,
                "open_gap_penalty": -1.0,
                "extend_gap_penalty": -0.5,
                "mode": "global"
            }
        }
    )

class PairwiseAlignmentResponse(BaseModel):
    score: float
    aligned_sequence1: str
    aligned_sequence2: str
    full_alignment_str: str
    parameters_used: Dict

class DownloadsMixin:
      
      def _download_routes_config(self):
        """Configure API routes, including Entrez downloads and Biopython tools."""

        @self.post(
            "/download/entrez",
            tags=["downloads"],
            summary="Download data from NCBI Entrez",
            operation_id="download_entrez_data",
            description="""
            Downloads data records from specified NCBI Entrez databases using Bio.Entrez.
            This endpoint is designed to be called by automated agents (like LLMs) or other services.

            **Critical Configuration:**
            The server hosting this API *must* have the `ENTREZ_EMAIL` environment variable set
            to a valid email address. NCBI requires this for Entrez queries to monitor usage
            and prevent abuse. Without it, NCBI may block requests.

            **Request Body Parameters:**
            - `ids` (List[str], required): A list of unique identifiers for the records to fetch
              from the specified Entrez database.
              Example: `["NM_000546.6", "AY123456.1"]`
            - `db` (DB_LITERAL, required): The target NCBI Entrez database.
              Common choices for sequences: 'nucleotide', 'protein'.
              Other examples: 'gene', 'pubmed', 'taxonomy'.
              For a comprehensive list of supported databases, refer to the `DB_LITERAL`
              type definition in the API schema (includes options like 'pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide', etc.).
              Ensure the `ids` provided are appropriate for the selected `db`.
            - `reftype` (Literal["fasta", "gb"], required): The desired format for the
              downloaded data.
                - "fasta": Returns data in FASTA format.
                - "gb": Returns data in GenBank format.
              Ensure the chosen `reftype` is compatible with the selected `db`.

            **Response:**
            - On success: Returns the downloaded data as a single raw string with a
              `Content-Type` of `text/plain`. The content directly corresponds to the
              data fetched from Entrez in the specified `reftype`.
            - On failure:
                - If NCBI Entrez returns an error (e.g., invalid ID, unsupported `db`/`reftype`
                  combination, rate limiting), an HTTPException with the corresponding
                  status code (e.g., 400, 404, 503) and details from NCBI will be raised.
                - For other unexpected server-side errors during the process, an
                  HTTPException with status code 500 will be raised.

            **Example LLM Usage:**
            An LLM agent intending to fetch the FASTA sequence for human TP53 mRNA (NM_000546.6)
            would construct a POST request to this endpoint with the following JSON body:
            ```json
            {{
                "ids": ["NM_000546.6"],
                "db": "nucleotide",
                "reftype": "fasta"
            }}
            ```
            The agent should then be prepared to handle a plain text response containing
            the FASTA sequence or an error object if the request fails.
            """,
        )
        def entrez_download(request: EntrezDownloadRequest):
            """
            Handles Entrez download requests.
            Uses the globally defined `get_entred` function.
            """
            try:
                downloaded_content = get_entred(
                    ids=request.ids,
                    db=request.db,
                    reftype=request.reftype
                )
                return downloaded_content
            except HTTPError as he:
                raise HTTPException(status_code=he.code, detail=f"NCBI Entrez Error ({he.code}): {he.reason}")
            except Exception as e:
                raise HTTPException(status_code=500, detail=f"An unexpected error occurred during Entrez download: {str(e)}")

        @self.post(
            "/tools/align/pairwise",
            response_model=PairwiseAlignmentResponse,
            tags=["biotools"],
            summary="Perform pairwise sequence alignment using Biopython",
            operation_id="perform_pairwise_alignment",
            description="""
            Performs a pairwise sequence alignment (global or local) using Biopython's PairwiseAligner.
            You can specify sequences and alignment scoring parameters.
            """
        )
        def pairwise_alignment_route(request: PairwiseAlignmentRequest):
            try:
                response = run_pairwise_alignment(request)
                return response
            except ValueError as ve:
                raise HTTPException(status_code=400, detail=str(ve))
            except Exception as e:
                raise HTTPException(status_code=500, detail=f"An unexpected error occurred during alignment: {str(e)}")
    
    

def run_pairwise_alignment(request: PairwiseAlignmentRequest) -> PairwiseAlignmentResponse:
    aligner = PairwiseAligner()
    aligner.match_score = request.match_score
    aligner.mismatch_score = request.mismatch_penalty
    aligner.open_gap_score = request.open_gap_penalty
    aligner.extend_gap_score = request.extend_gap_penalty
    aligner.mode = request.mode

    seq1 = Seq(request.sequence1)
    seq2 = Seq(request.sequence2)

    alignments = list(aligner.align(seq1, seq2))

    if not alignments:
        raise ValueError("No alignment could be produced with the given sequences and parameters.")

    best_alignment = alignments[0]
    
    aligned_seq1_str = str(best_alignment[0])
    aligned_seq2_str = str(best_alignment[1])

    return PairwiseAlignmentResponse(
        score=best_alignment.score,
        aligned_sequence1=aligned_seq1_str,
        aligned_sequence2=aligned_seq2_str,
        full_alignment_str=str(best_alignment),
        parameters_used=request.dict()
    )

def get_entred(ids: List[str], db: DB_LITERAL, reftype: Literal["fasta", "gb"]) -> str:
    """
    Downloads data from NCBI Entrez databases.

    This function uses Bio.Entrez to fetch data based on a list of IDs from a specified database.
    It requires an email address for NCBI, which should be set via the ENTREZ_EMAIL environment
    variable or it defaults to "default_email@example.com" if the variable is not set.

    Args:
        ids: A list of strings, where each string is an ID for the records to be downloaded
             from the specified Entrez database.
        db: The Entrez database to query. While many NCBI databases exist, this function's
            type hints suggest primary support for "protein" and "nucleotide".
            Example list of Entrez databases:
            ['pubmed', 'protein', 'nuccore', 'ipg', 'nucleotide', 'structure', 'genome',
            'annotinfo', 'assembly', 'bioproject', 'biosample', 'blastdbinfo', 'books',
            'cdd', 'clinvar', 'gap', 'gapplus', 'grasp', 'dbvar', 'gene', 'gds',
            'geoprofiles', 'medgen', 'mesh', 'nlmcatalog', 'omim', 'orgtrack', 'pmc',
            'proteinclusters', 'pcassay', 'protfam', 'pccompound', 'pcsubstance',
            'seqannot', 'snp', 'sra', 'taxonomy', 'biocollections', 'gtr'].
        reftype: The format of the returned data, e.g., 'fasta' or 'gb' (GenBank).
                 This function's type hints suggest primary support for "fasta" and "gb".

    Returns:
        A string containing the raw data downloaded from Entrez. For text-based formats
        like FASTA or GenBank, this is typically a string. For binary formats, it might be bytes
        (though this function is typed to return str).

    Raises:
        Various exceptions can be raised by Bio.Entrez.efetch if the request fails,
        such as HTTPError (e.g., for invalid IDs or incorrect db/reftype combinations),
        URLError (network issues), or other NCBI-specific errors.

    Example:
        >>> from pathlib import Path
        >>> import os
        >>> # CRITICAL: Set your email for NCBI Entrez access
        >>> # For example, in your shell: export ENTREZ_EMAIL="your.email@example.com"
        >>> # Fallback for example if not set (NCBI might block extensive use without email):
        >>> if not os.getenv("ENTREZ_EMAIL"):
        ...     print("Example: ENTREZ_EMAIL not set, using placeholder. NCBI usage may be restricted.")
        ...     os.environ["ENTREZ_EMAIL"] = "default_test@example.com" # For doctest only
        >>>
        >>> download_dir_example = Path("./entrez_example_output")
        >>> download_dir_example.mkdir(exist_ok=True)
        >>>
        >>> # Example: Fetching a nucleotide sequence in FASTA format
        >>> nucleotide_ids_ex = ["NM_000546.6"] # Human TP53 mRNA
        >>> if os.getenv("ENTREZ_EMAIL") != "default_test@example.com" or True: # Run if email is set
        ...     nucleotide_content = entrez_download(
        ...         where=download_dir_example,
        ...         ids=nucleotide_ids_ex,
        ...         db="nucleotide", # Defaulting to nucleotide
        ...         reftype="fasta"
        ...     )
        ...     print(f"Fetched nucleotide (first 60 chars): {nucleotide_content[:60].replace('\n', ' ')}...")
        >>>
        >>> # Example: Fetching a protein sequence (conceptual, check reftype for protein)
        >>> # To demonstrate a different database, let's use protein with an appropriate ID
        >>> protein_ids_ex = ["NP_000537.3"] # Human p53 protein
        >>> if os.getenv("ENTREZ_EMAIL") != "default_test@example.com" or True: # Run if email is set
        ...     protein_content_fasta = entrez_download(
        ...         where=download_dir_example,
        ...         ids=protein_ids_ex,
        ...         db="protein", 
        ...         reftype="fasta"
        ...     )
        ...     print(f"Fetched protein FASTA (first 60 chars): {protein_content_fasta[:60].replace('\n', ' ')}...")
    """
    Entrez.email = os.getenv("ENTREZ_EMAIL", "default_email@example.com")
    # Use the db and reftype parameters passed to the function
    handle = Entrez.efetch(db=db, id=ids, rettype=reftype)
    return handle.read()

if __name__ == "__main__":
    import os
    import time
    from Bio import Entrez
    from urllib.error import HTTPError # Import for specific error handling

    # --- NCBI Entrez Email Setup ---
    entrez_email_address = os.getenv("ENTREZ_EMAIL")
    if not entrez_email_address:
        print("WARNING: The ENTREZ_EMAIL environment variable is not set.")
        print("NCBI requires a valid email address for Entrez queries to prevent abuse.")
        print('Please set this variable, e.g., in your shell: `export ENTREZ_EMAIL="your.name@example.com"`')
        print("Using a default placeholder email for this run, which may lead to NCBI blocking requests.")
        Entrez.email = "default_cli_placeholder@example.com"
    else:
        Entrez.email = entrez_email_address
        print(f"Using Entrez email: {Entrez.email}")

    print("\n--- Querying NCBI Entrez for Database Information ---")
    print("Attempting to fetch the list of all available Entrez databases...")

    all_databases = []
    try:
        # Fetch the global list of databases
        handle = Entrez.einfo()
        record = Entrez.read(handle)
        handle.close()
        all_databases = record['DbList']
        print(f"Successfully fetched {len(all_databases)} database names.")
        # print(f"Databases found: {all_databases}") # Uncomment to see the full list immediately
    except Exception as e:
        print(f"FATAL ERROR: Could not fetch the list of Entrez databases: {type(e).__name__} - {e}")
        print("Cannot proceed without the database list.")
        exit(1) # Exit if we can't get the basic list

    print("\n--- Retrieving Rettpye/Retmode for each Database ---")
    print("This may take some time as it queries NCBI for each database...")

    db_info_results = {}
    databases_with_errors = []

    # Limit requests per second to NCBI
    requests_per_second_limit = 3 # NCBI recommends no more than 3 requests/sec without an API key
    delay_between_requests = 1.0 / requests_per_second_limit

    for db_name in all_databases:
        print(f"Querying info for database: '{db_name}'...")
        retries = 3 # Allow a few retries for transient network issues
        success = False
        for attempt in range(retries):
            try:
                handle = Entrez.einfo(db=db_name)
                record = Entrez.read(handle)
                handle.close()

                # Extract rettype and retmode lists
                db_info = record.get('DbInfo', {})
                rettype_list = db_info.get('RetTypeList', [])
                retmode_list = db_info.get('RetModeList', [])

                db_info_results[db_name] = {
                    'rettypes': rettype_list,
                    'retmodes': retmode_list
                }
                print(f"  -> Success: Found {len(rettype_list)} rettypes, {len(retmode_list)} retmodes.")
                success = True
                break # Exit retry loop on success

            except HTTPError as e:
                # Specific handling for HTTP errors (e.g., 4xx, 5xx) which might indicate db issues
                print(f"  -> HTTP Error for '{db_name}' (Attempt {attempt+1}/{retries}): {e.code} {e.reason}")
                # Don't retry immediately on persistent client/server errors like 400 or 404
                if e.code >= 400 and e.code < 500:
                     databases_with_errors.append((db_name, f"HTTP {e.code} {e.reason}"))
                     break
                # For server errors (5xx) or others, wait and retry
                time.sleep(delay_between_requests * (attempt + 1)) # Exponential backoff basic

            except Exception as e:
                # General error handling
                print(f"  -> ERROR querying '{db_name}' (Attempt {attempt+1}/{retries}): {type(e).__name__} - {e}")
                time.sleep(delay_between_requests * (attempt + 1)) # Basic backoff

        if not success and (db_name, f"Max retries reached or non-retryable error") not in databases_with_errors:
             databases_with_errors.append((db_name, "Failed after multiple attempts or non-retryable HTTP error"))
             print(f"  -> FAILED to get info for '{db_name}' after {retries} attempts.")

        # Pause between requests to respect NCBI rate limits
        time.sleep(delay_between_requests)

    print("\n--- Summary of Supported Rettpes and Retmodes per Database ---")

    # Sort databases alphabetically for consistent output
    sorted_db_names = sorted(db_info_results.keys())

    for db_name in sorted_db_names:
        info = db_info_results[db_name]
        print(f"\nDatabase: {db_name}")
        if info['rettypes']:
            print(f"  Supported RetTypes: {info['rettypes']}")
        else:
            print("  Supported RetTypes: None found or extraction failed.")
        if info['retmodes']:
            print(f"  Supported RetModes: {info['retmodes']}")
        else:
            print("  Supported RetModes: None found or extraction failed.")

    if databases_with_errors:
        print("\n--- Databases with Errors During Query ---")
        for db_name, error_msg in databases_with_errors:
            print(f"  Database: {db_name} - Error: {error_msg}")

    print("\nQuerying complete.")

