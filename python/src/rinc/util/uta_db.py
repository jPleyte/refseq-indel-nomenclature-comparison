'''
List the schemas in the UTA database.  
if the UTA_DB_URL environment variable is not this class connects to the remote uta.biocommons.org UTA database. 

Created on Jan 3, 2026

@author: pleyte
'''
import logging
import sys
from psycopg2.extras import DictCursor
import hgvs.dataproviders.uta as uta

# Setup logging to console
root = logging.getLogger()
root.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
root.addHandler(handler)
logger = logging.getLogger(__name__)

class UtaDb(object):
    '''
    classdocs
    '''
    def __init__(self):
        '''
        Constructor
        '''
        self._logger = logging.getLogger(__name__)
        self._hdp = None
    
    def __enter__(self):
        self._hdp = uta.connect()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._hdp.close()
        
    def get_latest_schema(self):
        return sorted(self.get_schemas(), reverse=True)[0]
    
    def get_schemas(self):
        """
        Print a list of the schema in the uta db
        """
        with self._hdp._conn.cursor(cursor_factory=DictCursor) as cursor:    
            cursor.execute("""SELECT schema_name 
                                FROM information_schema.schemata 
                               WHERE schema_name LIKE 'uta_%';
                            """)
            return [row[0] for row in cursor.fetchall()]
    
    def query(self, sql) -> list:
        """
        Execute query and return results
        """
        with self._hdp._conn.cursor(cursor_factory=DictCursor) as cursor:
            cursor.execute(sql)
            results = [ dict(x) for x in cursor ]
            return results

def main():
    with UtaDb() as uta_db:
        latest = uta_db.get_latest_schema()
        print(f"Available Schemas:") 
        for x in uta_db.get_schemas():
            print(f"  {x} {'(latest)' if x == latest else ''}") 

if __name__ == '__main__':
    main()