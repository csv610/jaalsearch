import streamlit as st
from tavily import TavilyClient
import os

# Ensure you have set the TAVILY_API_KEY environment variable
tavily_api_key = os.getenv("TAVILY_API_KEY")
if not tavily_api_key:
    st.error("Please set the TAVILY_API_KEY environment variable.")
    st.stop()

@st.cache_resource
def get_tavily_client():
    return TavilyClient(tavily_api_key)

def search_tavily(client, query, num_results=5):
    try:
        response = client.search(query, max_results=num_results)
        return response.get('results', [])
    except Exception as e:
        st.error(f"An error occurred while fetching search results: {e}")
        return []

st.set_page_config(page_title="Tavily AI-Powered Search", page_icon="🔍", layout="wide")

tavily_client = get_tavily_client()

st.sidebar.title("Tavily")
num_results = st.sidebar.slider("Number of results:", min_value=1, max_value=100, value=10)

search_query = st.text_input("Enter your search query:", on_change=lambda: st.session_state.update({'search_trigger': True}))

if search_query and (st.session_state.get('search_trigger') or st.button("Search")):
    st.session_state['search_trigger'] = False
    st.write(f"Searching for: {search_query}")
    results = search_tavily(tavily_client, search_query, num_results)
    
    if results:
        for i, result in enumerate(results, 1):
            st.subheader(f"Result {i}")
            st.write(f"**Link:** {result.get('url', 'N/A')}")
            st.write(f"**Summary:** {result.get('content', 'N/A')}")
            st.divider()
    else:
        st.write("No results found.")
