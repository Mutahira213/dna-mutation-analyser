import requests


class PDBFetcher:
    """
    Fetches 3D protein structure data from the RCSB Protein Data Bank (PDB).
    Steps:
      1. search_pdb()     — POST full-text search to find the PDB entry ID
      2. fetch_structure() — GET the .pdb file for 3D rendering with py3Dmol
    """

    SEARCH_URL   = "https://search.rcsb.org/rcsbsearch/v2/query"
    STRUCTURE_URL = "https://files.rcsb.org/download/{}.pdb"

    def __init__(self, protein_name: str):
        self.protein_name = protein_name
        self.pdb_id = None

    def search_pdb(self) -> str:
        """
        Searches PDB by protein name and returns the first matching PDB ID (e.g. '4HX3').
        Returns None if no results found. Raises ConnectionError on network failure.
        """
        query = {
            "query": {
                "type": "terminal",
                "service": "full_text",
                "parameters": {"value": self.protein_name}
            },
            "return_type": "entry",
            "request_options": {"paginate": {"start": 0, "rows": 1}}
        }
        try:
            response = requests.post(self.SEARCH_URL, json=query, timeout=10)
            response.raise_for_status()
            data = response.json()
            results = data.get("result_set", [])
            if results:
                self.pdb_id = results[0]["identifier"]
                return self.pdb_id
            return None
        except requests.RequestException as e:
            raise ConnectionError(f"PDB API error: {e}")

    def fetch_structure(self, pdb_id: str = None) -> str:
        """
        Downloads the full .pdb structure file as a string.
        Pass a pdb_id explicitly or rely on the one stored from search_pdb().
        Raises ConnectionError on network failure.
        """
        pid = pdb_id or self.pdb_id
        if not pid:
            raise ValueError("No PDB ID available. Run search_pdb() first.")
        try:
            url = self.STRUCTURE_URL.format(pid.upper())
            response = requests.get(url, timeout=15)
            response.raise_for_status()
            return response.text
        except requests.RequestException as e:
            raise ConnectionError(f"PDB structure fetch error: {e}")
