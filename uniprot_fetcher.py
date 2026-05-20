import requests


class UniProtFetcher:
    """
    Fetches protein annotation data from the UniProt REST API.
    Returns the protein's full name, biological function description,
    source organism, and amino acid sequence length.
    """

    SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

    def __init__(self, protein_name: str):
        self.protein_name = protein_name

    def search_protein(self) -> dict:
        """
        Queries UniProt for the protein name and returns the first raw result dict.
        Returns an empty dict if nothing found.
        Raises ConnectionError on network failure.
        """
        params = {
            "query": self.protein_name,
            "format": "json",
            "size": 1,
            "fields": "protein_name,organism_name,length,cc_function"
        }
        try:
            response = requests.get(self.SEARCH_URL, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()
            results = data.get("results", [])
            return results[0] if results else {}
        except requests.RequestException as e:
            raise ConnectionError(f"UniProt API error: {e}")

    def get_info(self) -> dict:
        """
        Parses the raw UniProt result and extracts:
          - name     : recommended or submitted protein name
          - function : biological function text from FUNCTION comment
          - organism : scientific name
          - length   : amino acid count
        Returns an empty dict if the protein is not found.
        """
        raw = self.search_protein()
        if not raw:
            return {}

        # --- Function description ---
        function = "Not available"
        for comment in raw.get("comments", []):
            if comment.get("commentType") == "FUNCTION":
                texts = comment.get("texts", [])
                if texts:
                    function = texts[0].get("value", "Not available")
                break

        # --- Organism ---
        organism = raw.get("organism", {}).get("scientificName", "Not available")

        # --- Sequence length ---
        length = raw.get("sequence", {}).get("length", "Not available")

        # --- Protein name (recommended > submitted fallback) ---
        name = "Not available"
        try:
            name = raw["proteinDescription"]["recommendedName"]["fullName"]["value"]
        except (KeyError, TypeError):
            try:
                name = raw["proteinDescription"]["submissionNames"][0]["fullName"]["value"]
            except (KeyError, TypeError):
                pass

        return {
            "name": name,
            "function": function,
            "organism": organism,
            "length": length
        }
