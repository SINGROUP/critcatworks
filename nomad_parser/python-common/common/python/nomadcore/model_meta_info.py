from __future__ import print_function
from past.builtins import cmp
from builtins import object
from sqlalchemy import Table, Column, Integer, String, Text, Unicode, UnicodeText, DateTime, ForeignKey, Boolean, SmallInteger, Enum, Float, and_
from sqlalchemy.orm import relationship, aliased
from sqlalchemy.orm.session import object_session
from sqlalchemy.sql.expression import insert, delete, select

import datetime, logging, json
from collections import namedtuple

from nomadcore.local_info_kinds import InfoKindEl, InfoKindEnv
from nomadcore.basic_meta_info import allMetaInfo, metaMetaInfo
from nomadcore.utils import goInteractive

from nomadcore.model_base import Base, useJson, compareSorted
if useJson:
    from nomadcore.model_base import JSON, JSONB

gidSize = 39
# hash  sizes(hexencoded) speed(with 64bit HW)
# md4    32    1 25%-50% more than md5
# md5    32    2 ~2 times sha256
# sha224 56    4
# sha256 64    4
# sha384 96    3 ~1.2-1.5 times sha224
# sha512 128   3
# sha1   40    not worth using, slower than sha256 in optimized implementations

# sha224 is the basis of the gid, then a character is added at the start depending on the
# a lower case character is prepended to indicate the method used to calculate the sha:
# u: uuid
# f: content of a file
# d: recursive sha of a directory
# p: a property sha
# i: dataTypeSha + standard formatted json data + propertyShas
# h: reproducible sha of hdf5 data
# j: sha of standard formatted json data
# t: timestamp

# ?: assimilation sha (prefix is the same as directory)
# s: a symbolic pseudo sha
# a: access group sha
# a: access control group sha
# T: timestampBlock
#
#
# A second character uppercase character might be prepended to build the egid:
# I: Info
# K: InfoKind
# A: AccessGroup or AccessControlGroup
# T: Trusted info (self generated, or vetted info)
# this lets the gid be trivially the sha of the content of the file without limitations
# or extra strings, while avoiding clashes arising from crafting files containing
# what is hashed in the other cases.
# The second character ensures that it is impossible to tamper with the separate
# access control no matter what one manags to store in the Info table.

info_kind_parents = Table('info_kind_parents', Base.metadata,
    Column('info_kind_gid', Integer, ForeignKey('info_kinds.gid'), nullable=False),
    Column('parent_gid', Integer, ForeignKey('info_kinds.gid'), nullable=False)
)
info_kind_parent = namedtuple('info_kind_parent', ['info_kind_gid', 'parent_gid'])

info_kind_anchestors = Table('info_kind_anchestors', Base.metadata,
    Column('info_kind_gid', Integer, ForeignKey('info_kinds.gid'), nullable=False),
    Column('anchestor_gid', Integer, ForeignKey('info_kinds.gid'), nullable=False)
)
info_kind_anchestor = namedtuple('info_kind_anchestor', ['info_kind_gid', 'anchestor_gid'])

class UnsetValue(object):
    pass

class InfoKind(Base):
    """Represents meta information, i.e. information on the information

    This is like a more powerful and structured set of tags.
    It is used both to define the type of information, and the pieces of information
    contained in it.

    It is always public, and versioning on them is also public (public_override)
    and thus has to be decided (vetted) by admins.

    There are 4 root types:
    - MetaType: implicit kind of the root InfoKind that describe InfoKind.kind.
        identified by a null InfoKind.kind
    - DocumentType: kind of InfoKinds describing the type of documents (Infos)
    - DocumentContentType: kind of InfoKinds describing the type of information contained
        in documents (Infos)
    - ConnectionType: kind of InfoKinds describing Connections
    - UnknownType: kind for the unknown InfoKinds kinds (only the name of them is relevant,
        to recover the correct kindStr)
    """
    __tablename__="info_kinds"
    gid = Column(String(gidSize), primary_key=True)
    name = Column(Unicode(128))
    description = Column(UnicodeText)
    kind_gid = Column(String(gidSize), ForeignKey('info_kinds.gid'), nullable=True)
    kind = relationship('InfoKind', foreign_keys=[kind_gid], remote_side = [gid])
    #meta_childs = relationship('InfoKind', foreign_keys=[kind_gid], remote_side=[gid])
    units = Column(Unicode(128))
    super_kinds = relationship("InfoKind",
                                                secondary=info_kind_parents,
                                                primaryjoin=gid==info_kind_parents.c.info_kind_gid,
                                                secondaryjoin=gid==info_kind_parents.c.parent_gid,
                                                backref="sub_kinds")
    super_kind_gid=Column(String(gidSize),ForeignKey('info_kinds.gid'), nullable=True)
    super_kind=relationship("InfoKind",foreign_keys=[super_kind_gid])
    cumulative_super_kinds = relationship("InfoKind",
                                                secondary=info_kind_anchestors,
                                                primaryjoin=gid==info_kind_anchestors.c.info_kind_gid,
                                                secondaryjoin=gid==info_kind_anchestors.c.anchestor_gid,
                                                backref="cumulative_sub_kinds")
    dtypeStr = Column(Unicode(128))
    repeats = Column(Boolean, nullable = True)
    if useJson:
        shape = Column(JSONB)
        extraArgs = Column(JSONB)
    else:
        shape_txt = Column(UnicodeText)
        extra_args_txt = Column(UnicodeText)
        def shape():
            doc = "The shape of array that stores this property (relevant for 'simple' properties)."
            def fget(self):
                if self.shape_txt:
                    return json.loads(self.shape_txt)
                return None
            def fset(self, value):
                self.shape_txt = json.dumps(value)
            return locals()
        shape = property(**shape())
        def extra_args():
            doc = "The extra_args of this InfoKinds."
            def fget(self):
                if self.extra_args_txt:
                    return json.loads(self.extra_args_txt)
                return None
            def fset(self, value):
                if value:
                    self.extra_args_txt = json.dumps(value).decode('utf-8')
            return locals()
        extra_args = property(**extra_args())
    public = Column(Boolean, default = True)
    public_override_gid = Column(String(gidSize), ForeignKey('info_kinds.gid'), nullable=True)
    public_override = relationship("InfoKind", foreign_keys=[public_override_gid])
    superseded_kinds = relationship("InfoKind", foreign_keys=[public_override_gid], remote_side=[gid])

    def super_gids():
        doc = "The super_gids property."
        def fget(self):
            res = [x.parent_gid for x in object_session(self).query(info_kind_parents.c.parent_gid).filter(info_kind_parents.c.info_kind_gid == self.gid).order_by(info_kind_parents.c.parent_gid).all()]
            if self.super_kind_gid:
                if self.super_kind_gid in res:
                    res.remove(self.super_kind_gid)
                res.insert(0, self.super_kind_gid)
            return res

        def fset(self, value):
            session = object_session(self)
            session.flush()
            oldValue = [i[0] for i in fget(self)]
            oldValue.sort()
            if value:
                self.super_kind_gid = value[0]
            else:
                self.super_kind_gid = None
            newValue = list(value)
            newValue.sort()
            toRemove, toAdd = compareSorted(oldValue, newValue)
            if toRemove:
                session.execute( delete(info_kind_parents).where(and_(info_kind_parents.c.gid == self.gid, authors_users.c.user_id.in_(toRemove))) )
            if toAdd:
                session.execute( insert(info_kind_parents).values([info_kind_parent(self.gid, gid) for gid in toAdd]) )
            session.commit()

        return locals()
    super_gids = property(**super_gids())

    def toDict(self, addExtraArgs = True, inlineExtraArgs = True , selfGid = True, subGids = True, precalculatedGid = False):
        res = {
            "description": self.description,
            "name": self.name,
            "super_names": self.super_names,
        }
        kString = self.kindStr
        if kString != "DocumentContentType":
            res["kindStr"] = kString
        if selfGid:
            res["gid"] = self.gid
        if subGids:
            res["super_gids"] = self.super_gids
        if self.units is not None:
            res["units"] = self.units
        if self.dtypeStr is not None:
            res["dtypeStr"] = self.dtypeStr
        if self.repeats is not None:
            res["repeats"] = self.repeats
        if self.shape is not None:
            res["shape"] = self.shape
        if addExtraArgs:
            if inlineExtraArgs:
                extraArgs = self.extra_args
                if extraArgs:
                    res.update(extraArgs)
            else:
                res["extra_args"] = self.extra_args
        return res

    def has_super_gid(self, gid):
        if self.super_kind_gid == gid:
            return True
        return object_session(self).query(info_kind_parents.c.parent_gid).filter(and_(info_kind_parents.c.info_kind_gid == self.gid, info_kind_parents.c.parent_gid == gid)).first() is not None

    def super_names():
        doc = "The super_names read only property."
        def fget(self):
            vals = object_session(self).query(InfoKind.gid, InfoKind.name).filter(InfoKind.gid == info_kind_parents.c.parent_gid).filter(info_kind_parents.c.info_kind_gid == self.gid).order_by(InfoKind.gid)
            res = []
            firstGid = self.super_kind_gid
            for val in vals:
                if firstGid == val.gid:
                    res.insert(0, val.name)
                else:
                    res.append(val.name)
            return res
        return locals()
    super_names = property(**super_names())

    def __init__(self, session = None, super_gids = None, super_names = None, kindStr = UnsetValue, **kwds):
        Base.__init__(self, **kwds)
        if session:
            session.add(self)
            session.flush()
        if super_gids:
            assert session, "session required for super_gids"
            self.super_gids = super_gids
        if kindStr is UnsetValue:
            if self.kind_gid == None:
                self.kind_gid = metaMetaInfo.gidOf('DocumentContentType')
        else:
            if self.kind_gid is None:
                self.kindStr=kindStr
            assert self.kindStr == kindStr, "kindStr is %s and but stored kind has name %s" % (kindStr, self.kindStr)
        if super_names is not None:
            if super_gids is None:
                if super_names:
                    raise Exception("super_gids required with super_names %s for %s" % (super_names, self.toDict()))
                else:
                    return
            if len(super_names) != len(super_gids):
                raise Exception("incompatible super_names and super_gids")
            # object_session(self).commit()
            sn=self.super_names
            if super_names:
                if super_names[0] != sn[0]:
                    raise Exception("first super type has different name (%r vs %r)" % (sn[0], super_names[0]))
                sn.sort()
                sn2 = list(super_names)
                sn2.sort()
                sn.sort()
                if sn != sn2:
                    raise Exception("unexpected type names in super_names (%r vs %r)" % (sn2, sn))

    _addedBasicTypes=False
    @classmethod
    def ensureBasicTypes(cls, session, force=False):
        cls._addedBasicTypes=True
        for metaInfo in allMetaInfo.values():
            cls.ensureInfoKindEnv(metaInfo, session)

    def kindStr():
        doc = "A string representing the kind of this InfoKind."
        def fget(self):
            if self.kind_gid is None:
                return "MetaType"
            else:
                ik = object_session(self).query(InfoKind.name).filter(InfoKind.gid == self.kind_gid).first()
                if ik:
                    return ik.name
                else:
                    return None

        def fset(self, value):
            if value == "MetaType":
                self.kind_gid = None
            else:
                if value in metaMetaInfo:
                    sha = metaMetaInfo.gidOf(value)
                    self.kind_gid = sha
                else:
                    knownIks = object_session(self).session.query(InfoKind).filter(InfoKind.name == value).filter(InfoKind.kind_gid == None).all()
                    if knownIks:
                        if len(knownIks)>1:
                            logging.getLogger("nomad.InfoKind").warn(
                                    "multiple %s InfoKinds in db detected %s when trying to add %s",
                                    ikEl.kindStr, [x.toInfoKindEl().toDict() for x in knownIks],ikEl)
                        knownIks.sort()
                        self.kind_gid = knownIks[0].gid
                    else:
                        logging.getLogger('nomad.InfoKind').warn('Unknown kind type %s cound not be resolved, replacing with an UnknownType',
                                ikEl.kindStr)
                        kindType = session.query(InfoKind).filter(and_(InfoKind.name == value, InfoKind.kind_gid == metaMetaInfo.gidOf('UnknownType'))).first()
                        if not kindType:
                            logging.getLogger('nomad.InfoKind').warn("will add unknown type for %s", ikEl.kindStr)
                            goInteractive(locals())
                            tmpEnv = InfoKindEnv()
                            tmpEnv.addInfoKindEl(metaMetaInfo.InfoKind('UnknownType'))
                            unkType = InfoKindEl(name = ikEl.kindStr, kindStr = 'UnknownType', description = 'auto inserted unknown type')
                            tmpEnv.addInfoKindEl(unkType)
                            self.kind = InfoKind(session = session, **unkType.toDict(env, inlineExtraArgs=False))
        return locals()
    kindStr = property(**kindStr())

    def egid():
            doc = "The egid property."
            def fget(self):
                    return "K" + self.gid
            def fset(self, value):
                    if len(value)!= 58 and len(value)!= 34:
                            raise Exception("Invalid length for egid value for InfoKind: {0!r}".format(value))
                    if value[0]!="K":
                        raise Exception("Invalid egid value for InfoKind: {0!r}".format(value))
                    self.gid = value[:gidSize]
            return locals()
    egid = property(**egid())

    @staticmethod
    def fromInfoKindEl(ikEl, env, session):
        gid = env.gidOf(ikEl.name)
        existingValue = session.query(InfoKind).get(gid)
        if existingValue is not None:
            return existingValue, False
        kindType = None
        if (ikEl.kindStr == '' or ikEl.kindStr == 'MetaType'):
            if not ikEl.name in metaMetaInfo:
                # ['DocumentType', 'DocumentContentType', 'AbstractDocumentContentType', 'ConnectionType', 'UnknownType']
                # allow???
                raise Exception('New root meta types not allowed, got {0}'.format(ikEl.toDict(env)))
        else:
            # kindStr should uniquely identify the type here we already handle the case that new root types are allowed
            basicKindTypeEl = metaMetaInfo.infoKindEl(ikEl.kindStr)
            envKindTypeEl = env.infoKindEl(ikEl.kindStr)
            if envKindTypeEl and not envKindTypeEl.kindStr in ["", "MetaType"]:
                envKindTypeEl = None
            if basicKindTypeEl and not basicKindTypeEl.kindStr in ["", "MetaType"]: # not a valid meta type, ignore
                basicKindTypeEl = None
            if basicKindTypeEl:
                kindType = session.query(InfoKind).get(metaMetaInfo.gidOf(basicKindTypeEl.name))
                if not kindType:
                    kindType = InfoKind(session = session, kind_gid = metaMetaInfo.gidOf(kindType.kindStr) if not kindType.kindStr in ["", "MetaType"] else None,
                            **basicKindTypeEl.toDict(metaMetaInfo, inlineExtraArgs = False))
                    session.commit()
            if envKindTypeEl is not None:
                if basicKindTypeEl:
                    # defined both in the env and as basic Type
                    if set(envKindTypeEl.super_names) != set(basicKindTypeEl.super_names):
                        logging.getLogger('nomad.InfoKind').warn('incompatible basic type %s vs %s, ignoring new type',
                            ikEl.toDict(), res.toDict())
                    elif gid != metaMetaInfo.gidOf(ikEl.kindStr):
                        logging.getLogger('nomad.InfoKind').warn('ignoring new kind type %s that has different gid but same name as a basic type %s',
                            kindType.toDict(), ikEl.toDict(env))
                else:
                    envKindType = session.query(InfoKind).get(env.gidOf(ikEl.kindStr))
                    if envKindType:
                        kindType = envKindType
                    else:
                        kindType = InfoKind(session = session, **envKindTypeEl.toDict(env, inlineExtraArgs = False))
            if not kindType:
                logging.getLogger('nomad.InfoKind').warn('Unknown kind type %s cound not be resolved, replacing with an UnknownType',
                            ikEl.kindStr)
                kindType = session.query(InfoKind).filter(InfoKind.name == ikEl.kindStr).filter(InfoKind.kind_gid == metaMetaInfo.gidOf('UnknownType')).first()
                if kindType is None:
                    print("will add unknown type for ", ikEl.kindStr)
                    goInteractive(locals())
                    tmpEnv = InfoKindEnv()
                    tmpEnv.addInfoKindEl(metaMetaInfo.infoKindEl('UnknownType'))
                    unkType = InfoKindEl(name = ikEl.kindStr, kindStr = 'UnknownType', description = 'auto inserted unknown type')
                    tmpEnv.addInfoKindEl(unkType)
                    assert metaMetaInfo.gidOf('UnknownType') == tmpEnv.gidOf('UnknownType')
                    kindType = InfoKind(session = session, kind_gid = tmpEnv.gidOf('UnknownType'), **unkType.toDict(tmpEnv, inlineExtraArgs=False))
                    session.commit()
        if ikEl.kindStr == 'DocumentType' or ikEl.kindStr == 'ConnectionType':
            # we would like the name to be enough to identify these, so we warn
            kindIK = aliased(InfoKind)
            dbIks = session.query(InfoKind).filter(InfoKind.name == ikEl.name).join(kindIK, InfoKind.kind).filter(kindIK.name == ikEl.kindStr).all()
            if len(dbIks)>0:
                knownIks=[]
                for dbIk in dbIks:
                    if not metaMetaInfo.gidOf('UnknownDocumentType') in dbIk.cumulative_super_kinds:
                        knownIks.append(dbIk)
                if knownIks:
                    knownIks.sort()
                    if len(knownIks)>1:
                        logging.getLogger("nomad.InfoKind").warn(
                                "multiple %s InfoKinds in db detected %s when trying to add %s",
                                ikEl.kindStr, [x.toInfoKindEl().toDict() for x in knownIks],ikEl)
                    res = knownIks[0]
                    superIK = aliased(InfoKind)
                    superNames = set(session.query(superIK.name).join(superIK, InfoKind.super_kinds).filter(InfoKind.gid == res.gid).all())
                    if superNames != set(ikEl.super_names):
                        logging.getLogger('nomad.InfoKind').warn('incompatible basic type %s vs %s',
                                ikEl.toDict(env), res.toDict())
                    else:
                        logging.getLogger('nomad.InfoKind').warn('found InfoKind %s when trying to add %s',
                                res.toDict(), ikEl.toDict(env))
        ikElDict = ikEl.toDict(env, inlineExtraArgs = False, selfGid = True, subGids = True)
        return InfoKind(session = session, kind = kindType, **ikElDict), True

    def toInfoKindEl(self):
        superKinds=self.super_kinds.all().order_by('name').values_list("name", flat=True)
        superKinds.sort()
        if self.super_kind:
            superKinds.remove(self.super_kind.name)
            superKinds.insert(0,self.super_kind.name)
        return InfoKindEl(name=self.name, description=self.description, kindStr=self.kindStr,
            units=self.units,super_names=superKinds)

    @classmethod
    def ensureInfoKindEnv(cls, env, session): 
        iks = env.sortedIKs()
        added = 0
        existing = 0
        for ik in iks:
            ki, wasAdded = cls.fromInfoKindEl(ik, env, session)
            if wasAdded:
                added += 1
            else:
                existing += 1
        return { 'added': added, 'existing': existing }

    def update_cumulative_super_kinds(self):
        newKinds = set()
        evaluatedKinds = set()
        toAdd=set([self])
        while len(toAdd)>0:
            k = toAdd.pop()
            if not k.publicOverride is None:
                k = k.publicOverride
            if k in evaluatedKinds:
                continue
            evaluatedKinds.add(k)
            newKinds.add(k)
            for k2 in k.super_kinds.all():
                if not k2.publicOverride is None:
                    k2 = k2.publicOverride
                if not k2 in evaluatedKinds:
                    toAdd.add(k2)
        toDel=set()
        for k in self.cumulative_super_kinds.all():
            if not k.publicOverride is None:
                k = k.publicOverride
            if k in newKinds:
                newKinds.remove(k)
            else:
                toDel.add(k)
        if len(toDel) == 0 and len(newKinds) == 0:
            return False
        for k in toDel:
            self.cumulative_super_kinds.remove(k)
        for k in newKinds:
            self.cumulative_super_kinds.add(k)
        return True

    def save(self, *argv, **kwds):
        if self.super_kind and not self.super_kind in self.super_kinds.all():
            self.super_kinds.add(self.super_kind)
        self.update_cumulative_super_kinds()
        super(InfoKind, self).save(*argv, **kwds)

    def __unicode__(self):
        return self.name

    def addPublicOverride(self, ik):
        if self.publicOverride==ik:
            return
        oldV=self.publicOverride
        self.publicOverride=ik
        self.save()
        for info in self.infos:
            info.update_has_info_of_kind()

    @staticmethod
    def connectionEgid():
        return "connectionTypeSK"

    @staticmethod
    def overridesConnectionEgid():
        return "overridesConnectionSK"

    @staticmethod
    def typeEgid():
        return "typeInfoKindSK"

    @staticmethod
    def propertyEgid():
        return "propertyTypeSK"

    def __cmp__(ik1, ik2):
        """orders first non overridden public, then overridden public, then by gid"""
        c = cmp(ik2.public, ik1.public)
        if c!= 0: return c
        c = cmp(ik2.public_override_gid is None, ik1.public_override_gid is None)
        if c!= 0: return c
        c = cmp(ik1.gid, ik2.gid)
        return c
