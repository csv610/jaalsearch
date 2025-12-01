# JAAlSearch - Multi-Engine Web Search Toolkit

A comprehensive, modular Python toolkit for searching across multiple web and scientific databases. Featuring a clean architecture with reusable search engine classes, CLI tools, and Streamlit web applications.

## Overview

JAAlSearch provides a unified interface for searching across multiple sources:
- **DuckDuckGo** - Privacy-focused general web search
- **Wikipedia** - Knowledge base search
- **PubMed** - Scientific and medical literature
- **SerpAPI** - Google search results (requires API key)
- **Tavily** - AI-powered search engine (requires API key)

## Project Structure

```
jaalsearch/
├── src/                          # Core module
│   ├── __init__.py
│   └── search_engines.py         # Reusable search classes
│
├── cli/                          # Command-line interfaces
│   ├── __init__.py
│   ├── duckgo.py                # DuckDuckGo CLI
│   ├── serp.py                  # SerpAPI CLI
│   └── tavily.py                # Tavily CLI
│
├── streamlit_apps/              # Web-based applications
│   ├── __init__.py
│   ├── duckgo.py                # DuckDuckGo web app
│   ├── serp.py                  # SerpAPI web app
│   ├── tavily.py                # Tavily web app
│   ├── wikipedia.py             # Wikipedia web app
│   └── pubmed.py                # PubMed web app
│
├── config/                      # Configuration
│   ├── __init__.py
│   └── settings.py              # App settings
│
├── utils/                       # Utilities and helpers
│   ├── __init__.py
│   └── helpers.py               # Helper functions
│
├── tests/                       # Test suite
│   ├── __init__.py
│   └── test_search_engines.py  # Search engine tests
│
├── docs/                        # Documentation
│   ├── REFACTORING_SUMMARY.md
│   ├── ARCHITECTURE.md
│   ├── USAGE_EXAMPLES.md
│   └── IMPLEMENTATION_SUMMARY.txt
│
├── requirements.txt             # Python dependencies
├── .gitignore                   # Git ignore rules
├── README.md                    # This file
└── REFACTORING_README.md        # Refactoring guide
```

## Quick Start

### Installation

```bash
git clone https://github.com/csv610/jaalsearch.git
cd jaalsearch
pip install -r requirements.txt
```

### Environment Variables (Optional)

```bash
export SERPAPI_API_KEY="your-key"
export TAVILY_API_KEY="your-key"
export PUBMED_EMAIL="your-email@example.com"
```

### CLI Usage

```bash
python -m cli.duckgo -q "python programming" -n 5
python -m cli.serp -q "machine learning" -n 10
python -m cli.tavily -q "artificial intelligence" -n 5
```

### Streamlit Apps

```bash
streamlit run streamlit_apps/duckgo.py
streamlit run streamlit_apps/serp.py
streamlit run streamlit_apps/tavily.py
streamlit run streamlit_apps/wikipedia.py
streamlit run streamlit_apps/pubmed.py
```

### Python Module

```python
from src import DuckDuckGoSearch, WikipediaSearch

ddg = DuckDuckGoSearch()
results = ddg.search("python tutorial", num_results=5)

wiki = WikipediaSearch()
results = wiki.search("Albert Einstein", num_results=3)
```

## Testing

```bash
python -m pytest tests/
python -m pytest tests/ --cov=src
```

## Documentation

See `docs/` folder for comprehensive documentation:
- REFACTORING_SUMMARY.md - Overview of refactoring
- ARCHITECTURE.md - Technical design
- USAGE_EXAMPLES.md - Usage examples
- IMPLEMENTATION_SUMMARY.txt - Metrics

## Features

✅ Reusable search engine classes  
✅ Multiple interfaces (CLI, Streamlit, Python API)  
✅ Clean architecture with separation of concerns  
✅ Comprehensive error handling and validation  
✅ Performance caching  
✅ Well-documented code and examples  
✅ Full test coverage  

## Environment Variables

| Variable | Purpose | Required |
|----------|---------|----------|
| `SERPAPI_API_KEY` | Google search | For SerpAPI |
| `TAVILY_API_KEY` | Tavily search | For Tavily |
| `PUBMED_EMAIL` | PubMed email | Optional |

## Dependencies

See `requirements.txt` for full list:
- requests
- duckduckgo-search
- wikipedia
- biopython
- tavily
- streamlit
- beautifulsoup4
- cachetools

## License

MIT License

---

Built with ❤️ for developers and researchers
