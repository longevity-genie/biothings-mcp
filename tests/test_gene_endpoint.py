import pytest
from fastapi.testclient import TestClient
from biothings_mcp.server import create_app
from pathlib import Path

@pytest.fixture
def client():
    """Fixture providing a FastAPI test client.
    
    This fixture creates a FastAPI test client that can be used to make HTTP requests
    to the API endpoints. The client is created using the create_app() function which
    initializes the FastAPI application with all the configured routes.
    """
    # Configure logging for tests
    project_root = Path(__file__).resolve().parents[1]  # Project root is one level up from tests dir
    log_dir = project_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)  # Ensure the directory exists

    app = create_app()
    return TestClient(app)

def test_get_gene_endpoint(client):
    """Test the /gene/{gene_id} endpoint.
    
    This test verifies that the endpoint correctly retrieves gene information
    for a given gene ID. It checks:
    1. The response status code is 200 (success)
    2. The response contains the correct gene ID
    3. The response contains the correct gene symbol
    4. The response contains the correct gene name
    
    The test uses the CDK2 gene (ID: 1017) as an example.
    """
    response = client.post("/gene/1017", json={})
    assert response.status_code == 200
    data = response.json()
    assert data["id"] == "1017"
    assert data["symbol"] == "CDK2"
    assert data["name"] == "cyclin dependent kinase 2"

def test_get_gene_with_fields_endpoint(client):
    """Test the /gene/{gene_id} endpoint with specific fields.
    
    This test verifies that the endpoint correctly filters the response to include
    only the requested fields. It checks:
    1. The response status code is 200 (success)
    2. The response contains only the requested fields (symbol and name)
    3. The response contains the required field (id)
    4. The field values are correct
    
    The test uses the CDK2 gene (ID: 1017) as an example and requests only
    the symbol and name fields.
    """
    response = client.post("/gene/1017", json={"fields": "symbol,name"})
    assert response.status_code == 200
    data = response.json()
    # Check that we have the requested fields
    assert "symbol" in data
    assert "name" in data
    assert data["symbol"] == "CDK2"
    assert data["name"] == "cyclin dependent kinase 2"
    # Check that we have the required fields
    assert "id" in data
    assert data["id"] == "1017"

def test_get_genes_endpoint(client):
    """Test the /genes endpoint for multiple genes.
    
    This test verifies that the endpoint correctly retrieves information for
    multiple genes in a single request. It checks:
    1. The response status code is 200 (success)
    2. The response contains the correct number of genes
    3. The response contains the correct gene IDs in the correct order
    
    The test uses two genes (CDK2 and CDK3) as examples.
    """
    # Make sure JSON body has the right structure
    response = client.post("/genes", json={
        "gene_ids": "1017,1018",
        "fields": None,  # Add default/expected fields
        "as_dataframe": False,
        "species": None
    })
    assert response.status_code == 200
    data = response.json()
    assert len(data) == 2
    assert data[0]["id"] == "1017"
    assert data[1]["id"] == "1018"

def test_query_genes_endpoint(client):
    """Test the /gene/query endpoint.
    
    This test verifies that the endpoint correctly performs a gene query and
    returns the expected results. It checks:
    1. The response status code is 200 (success)
    2. The response contains a hits array
    3. The hits array contains at least one result
    4. The first hit contains the correct gene information
    
    The test queries for the CDK2 gene using its symbol.
    """
    # Make sure JSON body has the right structure
    response = client.post("/gene/query", json={
        "q": "symbol:CDK2",
        "size": 1,
        "fields": None,
        "as_dataframe": False,
        "species": None,
        "sort": None,
        "skip": 0
    })
    assert response.status_code == 200
    data = response.json()
    assert "hits" in data
    assert len(data["hits"]) > 0
    hit = data["hits"][0]
    assert hit["symbol"] == "CDK2"
    assert hit["name"] == "cyclin dependent kinase 2"

def test_query_many_genes_endpoint(client):
    """Test the /gene/querymany endpoint.
    
    This test verifies that the endpoint correctly performs batch queries for
    multiple genes. It checks:
    1. The response status code is 200 (success)
    2. The response contains the correct number of results
    3. The response contains results for all queried genes
    
    The test queries for two genes (CDK2 and BRCA1) using their symbols.
    """
    # Make sure JSON body has the right structure
    response = client.post("/gene/querymany", json={
        "query_list": "CDK2,BRCA1",
        "scopes": "symbol",
        "size": 1,
        "fields": None,
        "as_dataframe": False,
        "species": None
    })
    assert response.status_code == 200
    data = response.json()
    assert len(data) == 2
    # Check that we got results for both queries
    symbols = {result["symbol"] for result in data}
    assert "CDK2" in symbols
    assert "BRCA1" in symbols

def test_metadata_endpoint(client):
    """Test the /gene/metadata endpoint.
    
    This test verifies that the endpoint correctly returns metadata about the
    gene database. It checks:
    1. The response status code is 200 (success)
    2. The response is a dictionary
    3. The response contains a stats field
    4. The stats field is a dictionary
    5. The stats field contains a total field
    6. The total field is an integer
    
    The metadata includes information about the total number of genes in the
    database and other statistics.
    """
    response = client.post("/gene/metadata", json={})
    assert response.status_code == 200
    data = response.json()
    assert isinstance(data, dict)
    assert "stats" in data
    assert isinstance(data["stats"], dict)
    assert "total" in data["stats"]
    assert isinstance(data["stats"]["total"], int)

def test_get_gene_ensembl_id(client):
    """Test the /gene/{gene_id} endpoint with an Ensembl ID.

    This test verifies that the endpoint correctly retrieves gene information
    using an Ensembl gene ID. It checks:
    1. The response status code is 200 (success)
    2. The response contains the correct gene ID (matches the queried Ensembl ID)
    3. Placeholder assertions for symbol and name (values need verification).

    The test uses the Ensembl ID ENSECAG00000002212 (likely from Cavia porcellus).
    """
    ensembl_id = "ENSECAG00000002212"
    response = client.post(f"/gene/{ensembl_id}", json={})
    assert response.status_code == 200
    data = response.json()
    # Note: The actual ID field returned by MyGene.info for Ensembl might be different
    # (e.g., it might resolve to an Entrez ID or keep the Ensembl ID in a specific field).
    # We now check if the queried ensembl_id exists within the 'ensembl.gene' field.
    assert "ensembl" in data
    assert "gene" in data["ensembl"]
    assert data["ensembl"]["gene"] == ensembl_id
    # TODO: Verify the expected symbol and name for this Ensembl ID and adjust assertions.
    # Check if standard fields are present, even if values aren't known
    assert "symbol" in data
    assert "name" in data
    assert "taxid" in data 