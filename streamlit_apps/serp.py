"""Google Search via SerpAPI Streamlit App"""

import streamlit as st
import sys
from pathlib import Path

# Add parent directory to path to import from src
sys.path.insert(0, str(Path(__file__).parent.parent))

from src import SerpAPISearch


@st.cache_resource
def get_search_engine():
    return SerpAPISearch()


st.set_page_config(layout="wide", page_title="Google Search", page_icon="🔍")

st.sidebar.title("SerpAPI")
num_results = st.sidebar.slider("Number of results:", min_value=1, max_value=100, value=10)

search_query = st.text_input("Enter your search query:")

if search_query.strip():
    st.write(f"🔍 Searching for: **{search_query}**")

    try:
        search_engine = get_search_engine()
        results = search_engine.search(search_query, num_results)

        if results:
            for i, result in enumerate(results, 1):
                st.subheader(f"Result {i}")
                st.write(f"**Title:** {result.get('title', 'No title')}")
                st.write(f"**Link:** {result.get('link', 'No link')}")
                st.write(f"**Summary:** {result.get('snippet', 'No description available')}")
                st.write("---")
        else:
            st.warning("No results found. Try a different search query.")
    except Exception as e:
        st.error(f"Error: {e}")
