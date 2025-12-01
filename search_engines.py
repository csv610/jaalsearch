"""
Reusable search engine classes for various search APIs.
Each class implements a consistent interface for performing searches.
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any
import os
import functools
import requests
from duckduckgo_search import DDGS
from tavily import TavilyClient
import wikipedia
import wikipedia.exceptions
from Bio import Entrez, Medline


class SearchEngine(ABC):
    """Abstract base class for all search engines."""

    @abstractmethod
    def search(self, query: str, num_results: int = 5) -> List[Dict[str, Any]]:
        """
        Perform a search query.

        Args:
            query: The search query string.
            num_results: Number of results to return.

        Returns:
            List of result dictionaries.
        """
        pass


class DuckDuckGoSearch(SearchEngine):
    """DuckDuckGo search implementation."""

    def __init__(self):
        self.ddgs = DDGS()

    def search(self, query: str, num_results: int = 5, safesearch: str = "Moderate") -> List[Dict[str, Any]]:
        """
        Perform a DuckDuckGo search.

        Args:
            query: The search query string.
            num_results: Number of results to return.
            safesearch: SafeSearch level ("Off", "Moderate", or "Strict").

        Returns:
            List of result dictionaries with 'title', 'href', and 'body' keys.
        """
        try:
            return list(self.ddgs.text(query, max_results=num_results, safesearch=safesearch))
        except Exception as e:
            raise Exception(f"DuckDuckGo search error: {e}")


class SerpAPISearch(SearchEngine):
    """Google search via SerpAPI implementation."""

    def __init__(self, api_key: str = None):
        self.api_key = api_key or os.getenv('SERPAPI_API_KEY')
        if not self.api_key:
            raise EnvironmentError("SERPAPI_API_KEY environment variable is not set.")

    def search(self, query: str, num_results: int = 5) -> List[Dict[str, Any]]:
        """
        Perform a Google search using SerpAPI.

        Args:
            query: The search query string.
            num_results: Number of results to return (1-100).

        Returns:
            List of result dictionaries with 'title', 'link', and 'snippet' keys.
        """
        if not query.strip():
            raise ValueError("Search query cannot be empty.")
        if not (1 <= num_results <= 100):
            raise ValueError("Number of results must be between 1 and 100.")

        try:
            params = {
                "q": query,
                "num": num_results,
                "api_key": self.api_key,
                "engine": "google"
            }
            response = requests.get('https://serpapi.com/search', params=params).json()
            return response.get('organic_results', [])
        except Exception as e:
            raise Exception(f"SerpAPI search error: {e}")


class TavilySearch(SearchEngine):
    """Tavily AI-powered search implementation."""

    def __init__(self, api_key: str = None):
        self.api_key = api_key or os.getenv("TAVILY_API_KEY")
        if not self.api_key:
            raise EnvironmentError("TAVILY_API_KEY environment variable is not set.")
        self.client = TavilyClient(self.api_key)

    def search(self, query: str, num_results: int = 5) -> List[Dict[str, Any]]:
        """
        Perform a Tavily search.

        Args:
            query: The search query string.
            num_results: Number of results to return.

        Returns:
            List of result dictionaries with 'url' and 'content' keys.
        """
        try:
            response = self.client.search(query, max_results=num_results)
            return response.get('results', [])
        except Exception as e:
            raise Exception(f"Tavily search error: {e}")


class WikipediaSearch(SearchEngine):
    """Wikipedia knowledge base search implementation."""

    def __init__(self):
        pass

    @functools.lru_cache(maxsize=128)
    def _cached_summary(self, title: str, sentences: int = 2) -> str:
        """Get cached summary for a Wikipedia page."""
        return wikipedia.summary(title, sentences=sentences)

    def search(self, query: str, num_results: int = 5) -> List[Dict[str, Any]]:
        """
        Perform a Wikipedia search.

        Args:
            query: The search query string.
            num_results: Number of results to return.

        Returns:
            List of result dictionaries with 'title', 'url', and 'summary' keys.
        """
        try:
            search_results = wikipedia.search(query, results=num_results)

            if not search_results:
                return []

            results = []
            for result in search_results:
                try:
                    page = wikipedia.page(result)
                    summary = self._cached_summary(result, sentences=2)

                    results.append({
                        "title": page.title,
                        "url": page.url,
                        "summary": summary
                    })
                except wikipedia.exceptions.DisambiguationError as e:
                    results.append({
                        "title": result,
                        "url": "",
                        "summary": f"Disambiguation page. Possible matches: {', '.join(e.options[:5])}"
                    })
                except wikipedia.exceptions.PageError:
                    results.append({
                        "title": result,
                        "url": "",
                        "summary": "Page not found"
                    })

            return results
        except Exception as e:
            raise Exception(f"Wikipedia search error: {e}")


class PubMedSearch(SearchEngine):
    """PubMed scientific literature search implementation."""

    def __init__(self, email: str = None):
        self.email = email or os.getenv("PUBMED_EMAIL", "csv610@gmail.com")
        Entrez.email = self.email

    def search(self, query: str, num_results: int = 5) -> List[Dict[str, Any]]:
        """
        Perform a PubMed search.

        Args:
            query: The search query string.
            num_results: Maximum number of results to return.

        Returns:
            List of result dictionaries with paper metadata including
            'title', 'authors', 'journal', 'pubdate', 'pmid', and 'abstract'.
        """
        try:
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=query, retmax=num_results)
            record = Entrez.read(handle)
            handle.close()

            # Fetch details for each result
            id_list = record["IdList"]
            if not id_list:
                return []

            handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
            records = Medline.parse(handle)

            results = []
            for record in records:
                results.append({
                    "title": record.get("TI", "No title available"),
                    "authors": ", ".join(record.get("AU", ["No authors listed"])),
                    "journal": record.get("JT", "No journal listed"),
                    "pubdate": record.get("DP", "No date available"),
                    "pmid": record.get("PMID", "No PMID available"),
                    "abstract": record.get("AB", "No abstract available")
                })

            return results
        except Exception as e:
            raise Exception(f"PubMed search error: {e}")
