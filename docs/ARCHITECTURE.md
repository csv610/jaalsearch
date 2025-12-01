# JAAlSearch Architecture

## Class Diagram

```
                    SearchEngine (ABC)
                    ├── search()

                           ▲
            ┌──────┬────────┼────────┬──────────┬──────────┐
            │      │        │        │          │          │
    DuckDuckGoSearch  SerpAPISearch  TavilySearch  WikipediaSearch  PubMedSearch
```

## Component Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                          Applications                             │
├──────────────────────┬──────────────────────┬──────────────────┤
│   CLI Applications   │  Streamlit Web Apps  │  Future: APIs    │
├──────────────────────┼──────────────────────┼──────────────────┤
│ cli_duckgo_search    │ sl_duckgo_search     │ (REST API)       │
│ cli_serp_search      │ sl_serp_search       │ (GraphQL API)    │
│ cli_tavily_search    │ sl_tavily_search     │                  │
│                      │ sl_wiki_search       │                  │
│                      │ sl_pubmed_search     │                  │
└──────────────────────┴──────────────────────┴──────────────────┘
                              ▲
                              │ imports
                              │ uses
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                   search_engines.py                              │
│  (Reusable Search Engine Classes)                               │
├─────────────────────────────────────────────────────────────────┤
│                                                                   │
│  SearchEngine (ABC)                                             │
│  ├── DuckDuckGoSearch()                                         │
│  ├── SerpAPISearch()                                            │
│  ├── TavilySearch()                                             │
│  ├── WikipediaSearch()                                          │
│  └── PubMedSearch()                                             │
│                                                                   │
└─────────────────────────────────────────────────────────────────┘
                              ▲
                              │ uses
                              │ integrates with
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                   External Services                              │
├──────────────────┬──────────────┬──────────────┬────────────────┤
│ DuckDuckGo API   │ SerpAPI      │ Tavily API   │ NCBI/PubMed   │ Wikipedia
│ (Free)           │ (Paid)       │ (Paid)       │ (Free)        │ (Free)
└──────────────────┴──────────────┴──────────────┴────────────────┘
```

## Data Flow

### CLI Example: DuckDuckGo Search

```
User Input (CLI args)
    │
    ├─ query: "python tutorial"
    ├─ num_results: 5
    └─ safesearch: "Moderate"

    ▼
cli_duckgo_search.py
    │
    ├─ Parse arguments
    ├─ Instantiate DuckDuckGoSearch()
    └─ Call search_engine.search()

    ▼
search_engines.py :: DuckDuckGoSearch
    │
    ├─ Validate inputs
    ├─ Call self.ddgs.text()
    └─ Return results

    ▼
Output (Formatted CLI results)
    ├─ Result 1: [title, link, summary]
    ├─ Result 2: [title, link, summary]
    └─ ...
```

### Streamlit Example: SerpAPI Search

```
User Input (Web UI)
    │
    ├─ Text input: "machine learning"
    └─ Slider: 10 results

    ▼
sl_serp_search.py
    │
    ├─ Get search engine (@st.cache_resource)
    ├─ Instantiate SerpAPISearch()
    └─ Call search_engine.search()

    ▼
search_engines.py :: SerpAPISearch
    │
    ├─ Validate query (not empty)
    ├─ Validate num_results (1-100)
    ├─ Make HTTPS request to https://serpapi.com/search
    └─ Return organic_results

    ▼
Output (Rendered Web UI)
    ├─ Result Card 1
    │  ├─ Title
    │  ├─ Link
    │  └─ Snippet
    ├─ Result Card 2
    └─ ...
```

## Class Interactions

### SearchEngine Base Class
```python
class SearchEngine(ABC):
    """Abstract base class for all search engines"""

    @abstractmethod
    def search(self, query: str, num_results: int = 5) -> List[Dict[str, Any]]:
        """Must be implemented by all subclasses"""
        pass
```

### Example: DuckDuckGoSearch
```python
class DuckDuckGoSearch(SearchEngine):
    def __init__(self):
        self.ddgs = DDGS()

    def search(self, query: str, num_results: int = 5, safesearch: str = "Moderate"):
        return list(self.ddgs.text(query, max_results=num_results, safesearch=safesearch))
```

### Example: SerpAPISearch
```python
class SerpAPISearch(SearchEngine):
    def __init__(self, api_key: str = None):
        self.api_key = api_key or os.getenv('SERPAPI_API_KEY')
        if not self.api_key:
            raise EnvironmentError("SERPAPI_API_KEY not set")

    def search(self, query: str, num_results: int = 5):
        # Validation
        if not query.strip():
            raise ValueError("Query cannot be empty")
        if not (1 <= num_results <= 100):
            raise ValueError("Results must be 1-100")

        # API request
        params = {
            "q": query,
            "num": num_results,
            "api_key": self.api_key,
            "engine": "google"
        }
        response = requests.get('https://serpapi.com/search', params=params).json()
        return response.get('organic_results', [])
```

## Result Format Consistency

Each search engine returns a list of dictionaries with engine-specific keys:

### DuckDuckGo
```python
[
    {
        "title": "...",
        "href": "...",
        "body": "..."
    }
]
```

### SerpAPI
```python
[
    {
        "title": "...",
        "link": "...",
        "snippet": "..."
    }
]
```

### Tavily
```python
[
    {
        "url": "...",
        "content": "..."
    }
]
```

### Wikipedia
```python
[
    {
        "title": "...",
        "url": "...",
        "summary": "..."
    }
]
```

### PubMed
```python
[
    {
        "title": "...",
        "authors": "...",
        "journal": "...",
        "pubdate": "...",
        "pmid": "...",
        "abstract": "..."
    }
]
```

## Dependency Injection & Configuration

### Environment Variables
```bash
SERPAPI_API_KEY=...         # Required for SerpAPI
TAVILY_API_KEY=...          # Required for Tavily
PUBMED_EMAIL=...            # Optional, defaults to csv610@gmail.com
```

### Optional Constructor Parameters
```python
search_engine = SerpAPISearch(api_key="custom-key")
search_engine = PubMedSearch(email="custom@email.com")
```

## Caching Strategy

### Streamlit Level (Application Cache)
```python
@st.cache_resource
def get_search_engine():
    return DuckDuckGoSearch()  # Reused across reruns
```

### Method Level (Function Cache)
```python
@functools.lru_cache(maxsize=128)
def _cached_summary(self, title: str, sentences: int = 2) -> str:
    return wikipedia.summary(title, sentences=sentences)
```

## Error Handling Flow

```
User Input
    │
    ▼
Validation (SearchEngine Class)
    │
    ├─ Invalid? ──► Raise ValueError
    │                   │
    │                   ▼
    │            CLI: print("Error: ...")
    │            Streamlit: st.error()
    │
    └─ Valid ─────► API Call
                        │
                        ├─ Success? ──► Return Results
                        │                   │
                        │                   ▼
                        │            Display to User
                        │
                        └─ Failed? ──► Raise Exception
                                           │
                                           ▼
                                    CLI: print("Error: ...")
                                    Streamlit: st.error()
```

## Extension Points

### Adding a New Search Engine

1. Create a new class in `search_engines.py`:
```python
class NewSearchEngine(SearchEngine):
    def __init__(self, api_key: str = None):
        self.api_key = api_key or os.getenv('NEW_API_KEY')
        self.client = NewClient(self.api_key)

    def search(self, query: str, num_results: int = 5):
        # Implementation here
        return results
```

2. Use it in CLI:
```python
from search_engines import NewSearchEngine
search_engine = NewSearchEngine()
results = search_engine.search(query, num_results)
```

3. Use it in Streamlit:
```python
from search_engines import NewSearchEngine

@st.cache_resource
def get_search_engine():
    return NewSearchEngine()

search_engine = get_search_engine()
results = search_engine.search(query, num_results)
```

No changes needed to existing files!

## Performance Considerations

### Current Optimizations
- **Caching**: `@st.cache_resource` prevents repeated instantiation
- **Summary Cache**: LRU cache for Wikipedia summaries
- **Lazy Loading**: Classes only imported when needed

### Future Optimizations
- Async/concurrent searches across multiple engines
- Connection pooling for API requests
- Result pagination for large result sets
- Incremental result loading in Streamlit
