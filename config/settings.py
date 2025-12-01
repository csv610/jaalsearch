"""
Configuration settings for JAAlSearch.

Environment variables and default settings.
"""

import os
from typing import Dict, Any

# Default number of search results
DEFAULT_RESULTS = 5
MAX_RESULTS = 100

# Supported search engines
SEARCH_ENGINES = {
    'duckduckgo': {
        'name': 'DuckDuckGo',
        'type': 'free',
        'requires_api_key': False,
        'privacy_level': 'high',
    },
    'wikipedia': {
        'name': 'Wikipedia',
        'type': 'free',
        'requires_api_key': False,
        'privacy_level': 'high',
    },
    'pubmed': {
        'name': 'PubMed',
        'type': 'free',
        'requires_api_key': False,
        'privacy_level': 'high',
    },
    'serpapi': {
        'name': 'SerpAPI (Google)',
        'type': 'paid',
        'requires_api_key': True,
        'privacy_level': 'medium',
        'env_var': 'SERPAPI_API_KEY',
    },
    'tavily': {
        'name': 'Tavily',
        'type': 'paid',
        'requires_api_key': True,
        'privacy_level': 'medium',
        'env_var': 'TAVILY_API_KEY',
    },
}

# API Keys from environment variables
API_KEYS: Dict[str, str] = {
    'serpapi': os.getenv('SERPAPI_API_KEY', ''),
    'tavily': os.getenv('TAVILY_API_KEY', ''),
}

# PubMed email configuration
PUBMED_EMAIL = os.getenv('PUBMED_EMAIL', 'csv610@gmail.com')

# SafeSearch levels for DuckDuckGo
SAFESEARCH_LEVELS = ['Off', 'Moderate', 'Strict']
DEFAULT_SAFESEARCH = 'Moderate'

# Streamlit configuration
STREAMLIT_CONFIG = {
    'layout': 'wide',
    'initial_sidebar_state': 'expanded',
}

# Logging configuration
LOGGING_CONFIG = {
    'level': os.getenv('LOG_LEVEL', 'INFO'),
    'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
}

# Development vs Production
DEBUG = os.getenv('DEBUG', 'False').lower() == 'true'
ENVIRONMENT = os.getenv('ENVIRONMENT', 'development')
