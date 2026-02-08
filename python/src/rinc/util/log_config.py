'''
Created on Jan 4, 2026

@author: pleyte
'''

class LogConfig(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        self.stdout_config = { 
            'version': 1,
            'disable_existing_loggers': False,
            'formatters': {
                'standard': { 
                    'format': '%(levelname)s: %(name)s::%(module)s:%(lineno)s: %(message)s'
                },
            },
            'handlers': {
                'default': {                     
                    'formatter': 'standard',
                    'class': 'logging.StreamHandler',
                    'stream': 'ext://sys.stdout'
                },
            },
            'loggers': { 
                '': {  # root logger
                    'level': 'WARNING',
                    'handlers': ['default'],
                    'propagate': False
                },
                'rinc': { 
                    'level': 'DEBUG',
                    'handlers': ['default'],
                    'propagate': False,
                },
                '__main__': { 
                    'level': 'DEBUG',
                    'handlers': ['default'],
                    'propagate': False,
                }
            }
        }
        
        self.file_config = { 
            'version': 1,
            'disable_existing_loggers': False,
            'formatters': {
                'standard': { 
                    'format': '%(asctime)s: %(levelname)s: %(name)s::%(module)s:%(lineno)s: %(message)s'
                },
            },
            'handlers': {
                'default': {                     
                    'formatter': 'standard',
                    'class': 'logging.FileHandler',
                    'filename': 'rinc.log', 
                    'mode': 'a'
                },
            },
            'loggers': { 
                '': {  # root logger
                    'level': 'WARNING',
                    'handlers': ['default'],
                    'propagate': False
                },
                'rinc': { 
                    'level': 'DEBUG',
                    'handlers': ['default'],
                    'propagate': False,
                },
                '__main__': { 
                    'level': 'DEBUG',
                    'handlers': ['default'],
                    'propagate': False,
                }
            }
        }