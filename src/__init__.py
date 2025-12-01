"""
JAAlSearch - Multi-Engine Web Search Toolkit

Core search engine classes and interfaces.
"""

from .search_engines import (
    SearchEngine,
    DuckDuckGoSearch,
    SerpAPISearch,
    TavilySearch,
    WikipediaSearch,
    PubMedSearch,
)

__version__ = "2.0.0"
__author__ = "JAAlSearch Contributors"

__all__ = [
    "SearchEngine",
    "DuckDuckGoSearch",
    "SerpAPISearch",
    "TavilySearch",
    "WikipediaSearch",
    "PubMedSearch",
]
