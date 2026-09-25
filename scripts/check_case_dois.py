#!/usr/bin/env python3
"""Check DOIs via Crossref API to verify they resolve correctly."""

import json
import urllib.request
import urllib.error
import sys

# DOIs from case meta.yaml files
dois = [
    # SMN2
    ("10.1128/MCB.26.4.1333-1346.2006", "Singh", "ISS-N1"),
    ("10.1056/NEJMoa1702752", "Finkel", "Nusinersen"),
    # CFTR
    ("10.1126/science.2475911", "Riordan", "CFTR gene"),
    ("10.1016/0888-7543(91)90503-7", "Zielenski", "CFTR structure"),
    # MAPT
    ("10.1074/jbc.274.21.15134", "Grover", "stem-loop"),
    ("10.1073/pnas.96.14.8229", "Varani", "NMR structure"),
    ("10.1038/31508", "Hutton", "FTDP-17"),
    # HCV
    ("10.1128/JVI.73.2.1165-1174.1999", "Honda", "domain II"),
    ("10.1038/nsb1004", "Lukavsky", "IRES structure"),
    # SARS-CoV-2
    ("10.1074/jbc.AC120.013449", "Kelly", "FSE mechanism"),
    ("10.1126/science.abf3546", "Bhatt", "FSE structure"),
]

print("Checking DOIs via Crossref API...")
print("=" * 80)

all_ok = True

for doi, expected_author, note in dois:
    url = f"https://api.crossref.org/works/{doi}"
    
    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            data = json.loads(response.read().decode())
        
        if 'message' in data:
            msg = data['message']
            title = msg.get('title', [''])[0] if 'title' in msg else 'NO TITLE'
            authors = msg.get('author', [])
            first_author = authors[0]['family'] if authors else 'NO AUTHOR'
            year = msg.get('published', {}).get('date-parts', [[None]])[0][0]
            
            author_match = expected_author.lower() in first_author.lower()
            
            status = "✓" if author_match else "✗"
            print(f"{status} {doi}")
            print(f"   Expected: {expected_author} ({note})")
            print(f"   Found: {first_author} et al. {year} — {title[:70]}...")
            print()
            
            if not author_match:
                all_ok = False
        else:
            print(f"✗ {doi} — No message in response")
            all_ok = False
    
    except urllib.error.HTTPError as e:
        print(f"✗ {doi} — HTTP {e.code}: {e.reason}")
        all_ok = False
    except Exception as e:
        print(f"✗ {doi} — Error: {e}")
        all_ok = False

print("=" * 80)

if all_ok:
    print("All DOIs verified successfully!")
    sys.exit(0)
else:
    print("Some DOIs failed verification")
    sys.exit(1)
