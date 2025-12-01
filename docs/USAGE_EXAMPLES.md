# JAAlSearch - Usage Examples

## Setup

### Environment Variables
```bash
export SERPAPI_API_KEY="your-serp-api-key"
export TAVILY_API_KEY="your-tavily-api-key"
# Optional: export PUBMED_EMAIL="your-email@example.com"
```

## CLI Usage

### DuckDuckGo Search

```bash
# Basic search with default settings (5 results, Moderate SafeSearch)
python cli_duckgo_search.py -q "python tutorial"

# Specify number of results
python cli_duckgo_search.py -q "machine learning" -n 10

# Disable SafeSearch filter
python cli_duckgo_search.py -q "hacking tutorials" -s "Off"

# Strict SafeSearch (recommended for all-ages content)
python cli_duckgo_search.py -q "educational content" -s "Strict"
```

### SerpAPI Search (Google Results)

```bash
# Basic Google search
python cli_serp_search.py -q "python programming"

# Get more results
python cli_serp_search.py -q "artificial intelligence" -n 20

# Search with custom number of results
python cli_serp_search.py -q "web development best practices" -n 15
```

### Tavily Search (AI-Powered)

```bash
# Basic Tavily search
python cli_tavily_search.py -q "latest technology trends"

# Tavily search with custom result count
python cli_tavily_search.py -q "machine learning frameworks" -n 10

# Search for specific topics
python cli_tavily_search.py -q "climate change research" -n 5
```

## Streamlit Web App Usage

### Run Individual Apps

```bash
# DuckDuckGo web interface
streamlit run sl_duckgo_search.py

# Google (SerpAPI) web interface
streamlit run sl_serp_search.py

# Tavily web interface
streamlit run sl_tavily_search.py

# Wikipedia web interface
streamlit run sl_wiki_search.py

# PubMed (scientific papers) web interface
streamlit run sl_pubmed_search.py
```

### Create a Multi-App Dashboard

Create `app.py`:
```python
import streamlit as st
import subprocess
import sys

st.set_page_config(page_title="JAAlSearch Dashboard", layout="wide")

st.title("🔍 JAAlSearch - Multi-Engine Search")

app_choice = st.sidebar.radio(
    "Select Search Engine",
    ["DuckDuckGo", "Google (SerpAPI)", "Tavily", "Wikipedia", "PubMed"]
)

if app_choice == "DuckDuckGo":
    st.write("### DuckDuckGo Search")
    exec(open("sl_duckgo_search.py").read())
elif app_choice == "Google (SerpAPI)":
    st.write("### Google Search")
    exec(open("sl_serp_search.py").read())
elif app_choice == "Tavily":
    st.write("### Tavily AI-Powered Search")
    exec(open("sl_tavily_search.py").read())
elif app_choice == "Wikipedia":
    st.write("### Wikipedia Search")
    exec(open("sl_wiki_search.py").read())
elif app_choice == "PubMed":
    st.write("### PubMed - Scientific Literature")
    exec(open("sl_pubmed_search.py").read())
```

Then run:
```bash
streamlit run app.py
```

## Programmatic Usage (Python)

### Using Search Engines in Your Code

#### DuckDuckGo
```python
from search_engines import DuckDuckGoSearch

# Initialize
ddg = DuckDuckGoSearch()

# Perform search
results = ddg.search("python tutorial", num_results=5, safesearch="Moderate")

# Process results
for result in results:
    print(f"Title: {result['title']}")
    print(f"URL: {result['href']}")
    print(f"Summary: {result['body']}")
    print("---")
```

#### SerpAPI (Google)
```python
from search_engines import SerpAPISearch

# Initialize (uses SERPAPI_API_KEY environment variable)
serp = SerpAPISearch()

# Or provide API key directly
serp = SerpAPISearch(api_key="your-api-key")

# Perform search
try:
    results = serp.search("machine learning", num_results=10)
    for result in results:
        print(f"Title: {result['title']}")
        print(f"URL: {result['link']}")
        print(f"Summary: {result['snippet']}")
except ValueError as e:
    print(f"Input error: {e}")
except Exception as e:
    print(f"Search error: {e}")
```

#### Tavily
```python
from search_engines import TavilySearch

# Initialize
tavily = TavilySearch()

# Perform search
results = tavily.search("latest AI breakthroughs", num_results=5)

for result in results:
    print(f"URL: {result['url']}")
    print(f"Content: {result['content']}")
    print("---")
```

#### Wikipedia
```python
from search_engines import WikipediaSearch

# Initialize
wiki = WikipediaSearch()

# Perform search
results = wiki.search("Albert Einstein", num_results=3)

for result in results:
    print(f"Title: {result['title']}")
    print(f"URL: {result['url']}")
    print(f"Summary: {result['summary']}")
    print("---")
```

#### PubMed
```python
from search_engines import PubMedSearch

# Initialize (uses PUBMED_EMAIL environment variable, defaults to csv610@gmail.com)
pubmed = PubMedSearch()

# Or provide email directly
pubmed = PubMedSearch(email="your-email@example.com")

# Perform search
try:
    results = pubmed.search("COVID-19 vaccines", num_results=5)

    for result in results:
        print(f"Title: {result['title']}")
        print(f"Authors: {result['authors']}")
        print(f"Journal: {result['journal']}")
        print(f"Date: {result['pubdate']}")
        print(f"PMID: {result['pmid']}")
        print(f"Abstract: {result['abstract'][:200]}...")  # First 200 chars
        print("---")
except Exception as e:
    print(f"Search error: {e}")
```

## Advanced Usage Patterns

### Batch Searching Multiple Queries

```python
from search_engines import DuckDuckGoSearch, SerpAPISearch, TavilySearch

queries = ["python", "javascript", "golang", "rust"]
ddg = DuckDuckGoSearch()

all_results = {}
for query in queries:
    try:
        results = ddg.search(query, num_results=3)
        all_results[query] = results
        print(f"✓ Found {len(results)} results for '{query}'")
    except Exception as e:
        print(f"✗ Error searching '{query}': {e}")
        all_results[query] = []

# Process results
for query, results in all_results.items():
    print(f"\n### Results for '{query}':")
    for result in results:
        print(f"  - {result['title']}")
```

### Multi-Engine Comparison

```python
from search_engines import DuckDuckGoSearch, WikipediaSearch

query = "artificial intelligence"

# Search with multiple engines
ddg = DuckDuckGoSearch()
wiki = WikipediaSearch()

ddg_results = ddg.search(query, num_results=3)
wiki_results = wiki.search(query, num_results=3)

print("=== DuckDuckGo Results ===")
for r in ddg_results:
    print(f"- {r['title']}: {r['href']}")

print("\n=== Wikipedia Results ===")
for r in wiki_results:
    print(f"- {r['title']}: {r['url']}")
```

### Error Handling

```python
from search_engines import SerpAPISearch

serp = SerpAPISearch()

try:
    # Empty query - will raise ValueError
    results = serp.search("", num_results=5)
except ValueError as e:
    print(f"Validation error: {e}")  # "Search query cannot be empty."

try:
    # Too many results - will raise ValueError
    results = serp.search("python", num_results=150)
except ValueError as e:
    print(f"Validation error: {e}")  # "Number of results must be between 1 and 100."

try:
    # API error - will raise Exception
    results = serp.search("test query", num_results=5)
except Exception as e:
    print(f"Search error: {e}")
```

### Filtering and Processing Results

```python
from search_engines import DuckDuckGoSearch

ddg = DuckDuckGoSearch()
results = ddg.search("python programming", num_results=20)

# Filter results by domain
python_org = [r for r in results if "python.org" in r.get('href', '')]
print(f"Found {len(python_org)} results from python.org")

# Filter by title keywords
tutorial_results = [r for r in results if 'tutorial' in r['title'].lower()]
print(f"Found {len(tutorial_results)} tutorial results")

# Extract and count keyword frequency
titles = [r['title'] for r in results]
keywords = {}
for title in titles:
    words = title.lower().split()
    for word in words:
        keywords[word] = keywords.get(word, 0) + 1

# Top 10 most common words
top_keywords = sorted(keywords.items(), key=lambda x: x[1], reverse=True)[:10]
print("Top keywords:", top_keywords)
```

### Combining Results from Multiple Sources

```python
from search_engines import DuckDuckGoSearch, WikipediaSearch

query = "climate change"

ddg = DuckDuckGoSearch()
wiki = WikipediaSearch()

# Collect results from both engines
all_results = []

try:
    ddg_results = ddg.search(query, num_results=3)
    all_results.extend(ddg_results)
except Exception as e:
    print(f"DuckDuckGo error: {e}")

try:
    wiki_results = wiki.search(query, num_results=3)
    all_results.extend(wiki_results)
except Exception as e:
    print(f"Wikipedia error: {e}")

# Remove duplicates and display
unique_titles = set()
for result in all_results:
    title = result.get('title', '')
    if title not in unique_titles:
        unique_titles.add(title)
        print(f"✓ {title}")

print(f"\nTotal unique results: {len(unique_titles)}")
```

### Saving Results to File

```python
import json
from search_engines import SerpAPISearch

serp = SerpAPISearch()
results = serp.search("machine learning", num_results=10)

# Save as JSON
with open("search_results.json", "w") as f:
    json.dump(results, f, indent=2)

# Save as CSV
import csv
with open("search_results.csv", "w", newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['title', 'link', 'snippet'])
    writer.writeheader()
    for result in results:
        writer.writerow(result)

print(f"Saved {len(results)} results")
```

## Testing & Validation

### Unit Test Example

```python
import unittest
from search_engines import DuckDuckGoSearch

class TestDuckDuckGoSearch(unittest.TestCase):
    def setUp(self):
        self.search = DuckDuckGoSearch()

    def test_search_returns_list(self):
        results = self.search.search("test", num_results=1)
        self.assertIsInstance(results, list)

    def test_search_with_high_safesearch(self):
        results = self.search.search("test", num_results=1, safesearch="Strict")
        self.assertIsInstance(results, list)

    def test_results_have_required_fields(self):
        results = self.search.search("python", num_results=1)
        if results:
            self.assertIn('title', results[0])
            self.assertIn('href', results[0])
            self.assertIn('body', results[0])

if __name__ == '__main__':
    unittest.main()
```

Run with:
```bash
python -m unittest test_search_engines.py
```
