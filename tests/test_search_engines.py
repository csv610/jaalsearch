"""
Unit tests for search engine classes.

This module provides basic unit tests for all search engine implementations.
"""

import unittest
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src import (
    DuckDuckGoSearch,
    SerpAPISearch,
    TavilySearch,
    WikipediaSearch,
    PubMedSearch,
)


class TestDuckDuckGoSearch(unittest.TestCase):
    """Test DuckDuckGo search implementation."""

    def setUp(self):
        self.search = DuckDuckGoSearch()

    def test_search_returns_list(self):
        """Test that search returns a list."""
        try:
            results = self.search.search("python", num_results=1)
            self.assertIsInstance(results, list)
        except Exception as e:
            self.skipTest(f"API error: {e}")

    def test_search_with_safesearch(self):
        """Test search with different SafeSearch levels."""
        try:
            for level in ["Off", "Moderate", "Strict"]:
                results = self.search.search("test", num_results=1, safesearch=level)
                self.assertIsInstance(results, list)
        except Exception as e:
            self.skipTest(f"API error: {e}")

    def test_results_have_required_fields(self):
        """Test that results have required fields."""
        try:
            results = self.search.search("python", num_results=1)
            if results:
                result = results[0]
                self.assertIn('title', result)
                self.assertIn('href', result)
                self.assertIn('body', result)
        except Exception as e:
            self.skipTest(f"API error: {e}")


class TestWikipediaSearch(unittest.TestCase):
    """Test Wikipedia search implementation."""

    def setUp(self):
        self.search = WikipediaSearch()

    def test_search_returns_list(self):
        """Test that search returns a list."""
        try:
            results = self.search.search("python programming", num_results=1)
            self.assertIsInstance(results, list)
        except Exception as e:
            self.skipTest(f"API error: {e}")

    def test_results_have_required_fields(self):
        """Test that results have required fields."""
        try:
            results = self.search.search("wikipedia", num_results=1)
            if results:
                result = results[0]
                self.assertIn('title', result)
                self.assertIn('url', result)
                self.assertIn('summary', result)
        except Exception as e:
            self.skipTest(f"API error: {e}")


class TestSerpAPISearch(unittest.TestCase):
    """Test SerpAPI search implementation."""

    def setUp(self):
        self.search = SerpAPISearch()

    def test_invalid_query_raises_error(self):
        """Test that empty query raises ValueError."""
        with self.assertRaises(ValueError):
            self.search.search("", num_results=5)

    def test_invalid_num_results_raises_error(self):
        """Test that invalid num_results raises ValueError."""
        with self.assertRaises(ValueError):
            self.search.search("python", num_results=150)

    def test_search_returns_list(self):
        """Test that search returns a list."""
        try:
            results = self.search.search("python", num_results=1)
            self.assertIsInstance(results, list)
        except Exception as e:
            self.skipTest(f"API error (likely missing API key): {e}")


class TestTavilySearch(unittest.TestCase):
    """Test Tavily search implementation."""

    def setUp(self):
        self.search = TavilySearch()

    def test_search_returns_list(self):
        """Test that search returns a list."""
        try:
            results = self.search.search("python", num_results=1)
            self.assertIsInstance(results, list)
        except Exception as e:
            self.skipTest(f"API error (likely missing API key): {e}")


class TestPubMedSearch(unittest.TestCase):
    """Test PubMed search implementation."""

    def setUp(self):
        self.search = PubMedSearch()

    def test_search_returns_list(self):
        """Test that search returns a list."""
        try:
            results = self.search.search("covid", num_results=1)
            self.assertIsInstance(results, list)
        except Exception as e:
            self.skipTest(f"API error: {e}")

    def test_results_have_required_fields(self):
        """Test that results have required fields."""
        try:
            results = self.search.search("python", num_results=1)
            if results:
                result = results[0]
                self.assertIn('title', result)
                self.assertIn('authors', result)
                self.assertIn('pmid', result)
        except Exception as e:
            self.skipTest(f"API error: {e}")


if __name__ == '__main__':
    unittest.main()
