from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import create_engine, exc
from sqlalchemy.orm import sessionmaker
import logging, sys

logger = logging.getLogger(__name__)
handler = logging.StreamHandler(stream=sys.stdout)
logger.setLevel(logging.INFO)
logger.addHandler(handler)

Base = declarative_base()
useNested = False
useJson = False

if useJson:
    from sqlalchemy.dialects.postgresql import JSON, JSONB

def compareSorted(oldValues, newValues):
    toAdd, toRemove = [], []
    iNew, iOld = 0, 0
    while iNew < len(newValues) and iOld < len(oldValues):
        if newValues[iNew] < oldValues[iOld]:
            toAdd.append(newValues[iNew])
            iNew += 1
        elif newValues[iNew] > oldValues[iOld]:
            toRemove.append(oldValues[iOld])
            iOld += 1
        else:
            iNew += 1
            iOld += 1
    toAdd += newValues[iNew:]
    toRemove += oldValues[iOld:]
    return (toRemove, toAdd)

def createEngine(engineStr='sqlite:///:memory:', echo = True):
    return create_engine(engineStr, echo = echo)

def createDB(engine):
    Base.metadata.create_all(engine)

def createSession(engine):
    return sessionmaker(bind=engine)()

def get_or_create(cls, session, defaults=None, **kwds):
    result = session.query(cls).filter_by(**kwds).first()
    if result:
        return result, True
    newVals=defaults
    if defaults is None:
        newVals={}
    newVals.update(kwds)
    result = cls(**newVals)
    if useNested:
        nestedSession = session.begin_nested()
        nestedSession.add(result)
        try:
            nestedSession.commit()
        except exc.IntegrityError:
            nestedSession.rollback()
            result = session.query(cls).filter_by(**kwds).one()
    else:
        session.add(result)
        session.flush()
    return result, True
