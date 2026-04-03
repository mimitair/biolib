class InterPro:
    """
    Class to interact with InterPro database
    """
    @staticmethod
    def getPDBAccessionsFromInterproIdentifier(interpro_identifier: str) -> list:
          """
          This function was copied from the interpro website
          https://www.ebi.ac.uk/interpro/result/download/#/structure/PDB/entry/InterPro/IPR029058/|accession
          Script was modified to take any IPR accession number, and to return the results as a list instead of writing to stdout.
          """
          # standard library modules
          import sys, errno, re, json, ssl
          from urllib import request
          from urllib.error import HTTPError
          from time import sleep
            
          print(f"Querying InterPro for PDB accessions associated with {interpro_identifier}.")  
          # MODIFIED: format the string to accept the given interpro accession
          BASE_URL =f"https://www.ebi.ac.uk:443/interpro/api/structure/PDB/entry/InterPro/{interpro_identifier}/?page_size=200"

          #disable SSL verification to avoid config issues
          context = ssl._create_unverified_context()

          next = BASE_URL
          last_page = False
          attempts = 0

          # MODIFIED: Empty list to store accession IDs
          result: list = []

          while next:
            try:
              req = request.Request(next, headers={"Accept": "application/json"})
              res = request.urlopen(req, context=context)
              # If the API times out due a long running query
              if res.status == 408:
                # wait just over a minute
                sleep(61)
                # then continue this loop with the same URL
                continue
              elif res.status == 204:
                #no data so leave loop
                break
              payload = json.loads(res.read().decode())
              next = payload["next"]
              attempts = 0
              if not next:
                last_page = True
            except HTTPError as e:
              if e.code == 408:
                sleep(61)
                continue
              else:
                # If there is a different HTTP error, it wil re-try 3 times before failing
                if attempts < 3:
                  attempts += 1
                  sleep(61)
                  continue
                else:
                  sys.stderr.write("LAST URL: " + next)
                  raise e

            # This is the part that returns the identifiers, notice that it is still in the While loop. 
            for i, item in enumerate(payload["results"]):
              # MODIFIED: Storing accessions in a list instead of writing to stdout:
              result.append(item["metadata"]["accession"])
              #sys.stdout.write(item["metadata"]["accession"] + "\n")
     
            # Don't overload the server, give it time before asking for more
            if next:
              sleep(1)

          # MODIFIED: return the result
          return result

