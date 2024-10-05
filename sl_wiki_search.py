import streamlit as st
import wikipedia
import functools
import cachetools.func
from bs4 import BeautifulSoup

@functools.lru_cache(maxsize=128)
def cached_summary(title, sentences=2):
    return wikipedia.summary(title, sentences=sentences)

def search_wikipedia(query, num_results=5):
    try:
        # Search for the query
        search_results = wikipedia.search(query, results=num_results)
        
        if not search_results:
            return []
        
        results = []
        for result in search_results:
            try:
                # Get a summary of the page
                page = wikipedia.page(result)
                summary = cached_summary(result, sentences=2)
                
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
    
    except wikipedia.exceptions.WikipediaException as e:
        st.error(f"A Wikipedia error occurred: {str(e)}")
    except ValueError as e:
        st.error(f"A value error occurred: {str(e)}")
    except Exception as e:
        st.error(f"An unexpected error occurred: {str(e)}")
    
    return []

# Streamlit app
st.sidebar.title("Wikipedia Search")

# Sidebar
num_results = st.sidebar.slider("Number of results", 1, 50, 10)

# Main area
query = st.text_input("Enter your search query and press Enter:")

if query:
    with st.spinner("Searching Wikipedia..."):
        results = search_wikipedia(query, num_results)
    
    if results:
        for result in results:
            st.markdown(f"### {result['title']}")
            if result['url']:
                st.markdown(f"[Link to Wikipedia page]({result['url']})")
            st.write(result['summary'])
            st.divider()
    else:
        st.warning(f"No results found for '{query}'.")
