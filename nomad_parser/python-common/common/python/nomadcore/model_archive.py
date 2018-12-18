from __future__ import print_function
from past.builtins import cmp
from builtins import object
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Table, Column, Integer, String, Text, Unicode, UnicodeText, DateTime, ForeignKey, Boolean, SmallInteger, Enum, Float
from sqlalchemy import create_engine, exc, Index, UniqueConstraint, and_
from sqlalchemy.orm import relationship, aliased
from sqlalchemy.orm.session import object_session
import datetime
import logging
import random
from nomadcore import storage_document, compact_sha 
from nomadcore.local_info_kinds import InfoKindEl, InfoKindEnv
from nomadcore.basic_info_kinds import basicInfoKinds
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from nomadcore.model_user import AccessGroup
from nomadcore.model_base import Base, useJson, compareSorted
from sqlalchemy.sql.expression import insert, delete, select
if useJson:
    from nomadcore.model_base import JSON, JSONB
from collections import namedtuple
import json

import readline # optional, will allow Up/Down/History in the console
import code

def goInteractive(locals):
    vars = globals().copy()
    vars.update(locals)
    shell = code.InteractiveConsole(vars)
    shell.interact()

Sha=compact_sha.sha224

def createEngine(engineStr='sqlite:///:memory:'):
    return create_engine(engineStr, echo=True)

def createDB(engine):
    Base.metadata.create_all(engine)

def createSession(engine):
    return sessionmaker(bind=engine)()

# hash  sizes(hexencoded) speed(with 64bit HW)
# md4    32    1 25%-50% more than md5
# md5    32    2 ~2 times sha256
# sha224 56    4
# sha256 64    4
# sha384 96    3 ~1.2-1.5 times sha224
# sha512 128   3
# sha1   40    not worth using, slower than sha256 in optimized implementations

# sha256 is the basis of the gid, then a character is added at the start depending on the
# method used to chalculate the sha:
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
# A second character is added to build the egid:
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
    gid = Column(String(57), primary_key=True)
    name = Column(Unicode(128))
    description = Column(UnicodeText)
    kind_gid = Column(String(57), ForeignKey('info_kinds.gid'), nullable=True)
    kind = relationship('InfoKind', foreign_keys=[kind_gid], remote_side = [gid])
    #meta_childs = relationship('InfoKind', foreign_keys=[kind_gid], remote_side=[gid])
    units = Column(Unicode(128))
    super_kinds = relationship("InfoKind",
                                                secondary=info_kind_parents,
                                                primaryjoin=gid==info_kind_parents.c.info_kind_gid,
                                                secondaryjoin=gid==info_kind_parents.c.parent_gid,
                                                backref="sub_kinds")
    super_kind_gid=Column(String(57),ForeignKey('info_kinds.gid'), nullable=True)
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
                    self.extra_args_txt = json.dumps(value)
            return locals()
        extra_args = property(**extra_args())
    public = Column(Boolean, default = True)
    public_override_gid = Column(String(57), ForeignKey('info_kinds.gid'), nullable=True)
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
                self.kind_gid = basicInfoKinds.gidOf('DocumentContentType')
        else:
            if self.kind_gid is None:
                self.kindStr=kindStr
            assert self.kindStr == kindStr, "kindStr is %s and but stored kind has name %s" % (kindStr, self.kindStr)
        if super_names is not None:
            if super_gids is None and super_names :
                raise Exception("super_gids required with super_names")
            if len(super_names) != len(super_gids):
                raise Exception("incompatible super_names and super_gids")
            # object_session(self).commit()
            sn=self.super_names
            if super_names:
                print("super_names:", super_names, "sn:", sn, "super_gids", super_gids)
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
        cls.ensureInfoKindEnv(basicInfoKinds, session)

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
                if value in basicInfoKinds:
                    sha = basicInfoKinds.gidOf(value)
                    self.kind_gid = sha
                else:
                    dbIks = object_session(self).session.query(InfoKind).filter(InfoKind.name == value).all()
                    knownIks=[]
                    if len(dbIks)>0:
                        for dbIk in dbIks:
                            if not basicInfoKinds.gidOf('UnknownDocumentType') in dbIk.cumulative_super_kinds:
                                    knownIks.append(dbIk)
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
                        kindType = session.query(InfoKind).filter(and_(InfoKind.name == value, InfoKind.kind_gid == basicInfoKinds.gidOf('UnknownType'))).first()
                        if not kindType:
                            print("will add unknon type for ", ikEl.kindStr)
                            goInteractive(locals())
                            tmpEnv = InfoKindEnv()
                            tmpEnv.addInfoKindEl(basicInfoKinds.InfoKind('UnknownType'))
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
                    self.gid = value[:57]
            return locals()
    egid = property(**egid())

    @staticmethod
    def fromInfoKindEl(ikEl, env, session):
        gid = env.gidOf(ikEl.name)
        existingValue = session.query(InfoKind).get(gid)
        if existingValue is not None:
            return existingValue
        kindType = None
        if (ikEl.kindStr == '' or ikEl.kindStr == 'MetaType'):
            if not ikEl.name in ['DocumentType', 'DocumentContentType', 'AbstractDocumentContentType', 'ConnectionType', 'UnknownType']:
                #allow???
                raise Exception('New root meta types not allowed, got {0}'.format(ikEl.toDict(env)))
        else:
            # kindStr should uniquely identify the type here we already handle the case that new root types are allowed
            basicKindTypeEl = basicInfoKinds.infoKindEl(ikEl.kindStr)
            envKindTypeEl = env.infoKindEl(ikEl.kindStr)
            if envKindTypeEl and not envKindTypeEl.kindStr in ["", "MetaType"]:
                envKindTypeEl = None
            if basicKindTypeEl and not basicKindTypeEl.kindStr in ["", "MetaType"]: # not a valid meta type, ignore
                basicKindTypeEl = None
            if basicKindTypeEl:
                kindType = session.query(InfoKind).get(basicInfoKinds.gidOf(basicKindTypeEl.name))
                if not kindType:
                    kindType = InfoKind(session = session, kind_gid = basicInfoKinds.gidOf(kindType.kindStr) if not kindType.kindStr in ["", "MetaType"] else None,
                            **basicKindTypeEl.toDict(basicInfoKinds, inlineExtraArgs = False))
                    session.commit()
            if envKindTypeEl is not None:
                if basicKindTypeEl:
                    # defined both in the env and as basic Type
                    if set(envKindTypeEl.super_names) != set(basicKindTypeEl.super_names):
                        logging.getLogger('nomad.InfoKind').warn('incompatible basic type %s vs %s, ignoring new type',
                            ikEl.toDict(), res.toDict())
                    elif gid != basicInfoKinds.gidOf(ikEl.kindStr):
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
                kindType = session.query(InfoKind).filter(InfoKind.name == ikEl.kindStr).filter(InfoKind.kind_gid == basicInfoKinds.gidOf('UnknownType')).first()
                if kindType is None:
                    print("will add unknon type for ", ikEl.kindStr)
                    goInteractive(locals())
                    tmpEnv = InfoKindEnv()
                    tmpEnv.addInfoKindEl(basicInfoKinds.infoKindEl('UnknownType'))
                    unkType = InfoKindEl(name = ikEl.kindStr, kindStr = 'UnknownType', description = 'auto inserted unknown type')
                    tmpEnv.addInfoKindEl(unkType)
                    assert basicInfoKinds.gidOf('UnknownType') == tmpEnv.gidOf('UnknownType')
                    kindType = InfoKind(session = session, kind_gid = tmpEnv.gidOf('UnknownType'), **unkType.toDict(tmpEnv, inlineExtraArgs=False))
                    session.commit()
        if ikEl.kindStr == 'DocumentType' or ikEl.kindStr == 'ConnectionType':
            # we would like the name to be enough to identify these, so we warn
            kindIK = aliased(InfoKind)
            dbIks = session.query(InfoKind).filter(InfoKind.name == ikEl.name).join(kindIK, InfoKind.kind).filter(kindIK.name == ikEl.kindStr).all()
            if len(dbIks)>0:
                knownIks=[]
                for dbIk in dbIks:
                    if not basicInfoKinds.gidOf('UnknownDocumentType') in dbIk.cumulative_super_kinds:
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
        ikElDict = ikEl.toDict(env, inlineExtraArgs = False)
        return InfoKind(session = session, kind=kindType, **ikElDict)

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
        for ik in iks:
            cls.fromInfoKindEl(ik, env, session)

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
#class InfoAnnotation(models.Model):
#  name = models.CharField(max_length=128)
#  description = models.TextField()
#  gid = models.CharField(max_length=57, db_index=True)
#  contains = models.ManyToManyField("self", related_name="contained_in", symmetrical=False)
#  related_info = models.ManyToManyField("self")
#  derived_from = models.ManyToManyField("self", related_name="derived_info", symmetrical=False)
#  to_index = models.TextField()
#  access_group = models.ForeignKey(AccessGroup, related_name="infos")
#  rejecters = models.ManyToManyField(settings.AUTH_USER_MODEL, related_name = "rejected_infos")

# class PublicInfoQuerySet(query.QuerySet):
#   'Query set with extra methods to select the various public visibility levels'
#   def visibleNoLogin(self):
#     "Info that is accessible without login"
#     return self.filter(access_group=AccessGroup.noLoginGroupId())
#   def activeNoLogin(self):
#     "Info that is active and visible without login"
#     return self.filter(active_group=AccessGroup.noLoginGroupId())
#   def visiblePublic(self):
#     "Info that is accessible to anybody after a login"
#     return self.filter(access_group__id__in=[AccessGroup.noLoginGroupId(),AccessGroup.publicGroupId()])
#   def activePublic(self):
#     "Info that is accessible to anybody after a login"
#     return self.filter(active_group__id__in=[AccessGroup.noLoginGroupId(),AccessGroup.publicGroupId()])
# 
# class InfoQuerySet(PublicInfoQuerySet):
#   def visiblePrivateForUser(self, user):
#     userGroups=AccessGroup.groupsForUser(user)
#     return self.filter(private_access_group__id__in=userGroups)
#   def activePrivateForUser(self, user):
#     userGroups=AccessGroup.groupsForUser(user)
#     return self.filter(private_active_group__id__in=userGroups).exclude(rejecters=user)
#   def visibleRestrictedorUser(self, user):
#     userGroups=AccessGroup.groupsForUser(user)
#     return self.filter(access_group__id__in=userGroups)
#   def activeRestrictedForUser(self, user):
#     userGroups=AccessGroup.groupsForUser(user)
#     return self.filter(active_group__id__in=userGroups).exclude(rejecters=user)
#   def visibleForUser(self, user):
#     userGroups=AccessGroup.groupsForUser(user)
#     userGroups.append(AccessGroup.publicGroupId())
#     userGroups.append(AccessGroup.noLoginGroupId())
#     return self.filter(access_group__id__in=userGroups)
#   def activeForUser(self, user):
#     userGroups=AccessGroup.groupsForUser(user)
#     userGroups.append(AccessGroup.publicGroupId())
#     userGroups.append(AccessGroup.noLoginGroupId())
#     return self.filter(active_group__id__in=userGroups).exclude(rejecters=user)
# 
# class RestrictedInfoQuerySet(PublicInfoQuerySet):
#   def __init__(self, model, accessGroups, rejecter=None, **kwds):
#     super(RestrictedInfoQuerySet,self).__init__(model,**kwds)
#     self.accessGroups=accessGroups
#     self.rejecter=rejecter
#   def activePrivate(self):
#     'informations that are active for one of the listed access groups'
#     res=self.filter(private_active_group__id__in=self.accessGroups)
#     if self.rejecter:
#       res=res.exclude(rejecters=self.rejecter)
#     return res
#   def visiblePrivate(self):
#     'informations that are visible for one of the listed access groups'
#     return self.filter(private_access_group__id__in=self.accessGroups)
#   def activeRestricted(self):
#     'all active informations that have been added to one of the listed access groups and is not public'
#     res=self.filter(active_group__id__in=self.accessGroups)
#     if self.rejecter:
#       res=res.exclude(rejecters=self.rejecter)
#     return res
#   def visibleRestricted(self):
#     'all visible informations that have been added to one of the listed access groups and is not public'
#     return self.filter(access_group__id__in=self.accessGroups)
#   def active(self):
#     accessGroups=set(self.accessGroups)
#     accessGroups.append(AccessGroup.publicGroupId())
#     accessGroups.append(AccessGroup.noLoginGroupId())
#     return self.filter(active_group__id__in=accessGroups)
#   def visible(self):
#     accessGroups=set(self.accessGroups)
#     accessGroups.append(AccessGroup.publicGroupId())
#     accessGroups.append(AccessGroup.noLoginGroupId())
#     return self.filter(access_group__id__in=accessGroups)
# 
# class InfoManager(models.Manager):
#   def get_queryset(self):
#     return InfoQuerySet(self.model, using=self._db)
# 
#   def forUser(self, user):
#     return RestrictedInfoQuerySet(self.model, AccessGroup.groupsForUser(user),
#         rejecter=user, using=self._db)
# 
# class AccessGroupsInfoManager(models.Manager):
#   def __init__(self, accessGroups, rejecter=None, **kwds):
#     super(UserInfoManager, self).__init__(**kwds)
#     self._accessGroups=accessGroups
#     self.rejecter=rejecter
#   def get_queryset(self):
#     return RestrictedInfoQuerySet(self.model, self.accessGroups,
#         rejecter=self.rejecter, using=self._db)
# 
# class UserInfoManager(models.Manager):
#   def __init__(self, user, **kwds):
#     super(UserInfoManager, self).__init__(**kwds)
#     self.user=user
#     self._accessGroups=None

info_base_info_kinds = Table('info_base_info_kinds', Base.metadata,
    Column('info_gid', Integer, ForeignKey('infos.gid'), nullable=False),
    Column('info_kind_gid', Integer, ForeignKey('info_kinds.gid'), nullable=False)
)

info_info_kinds = Table('info_info_kinds', Base.metadata,
    Column('info_gid', Integer, ForeignKey('infos.gid'), nullable=False),
    Column('info_kind_gid', Integer, ForeignKey('info_kinds.gid'), nullable=False)
)

info_rejecters = Table('info_rejecters', Base.metadata,
    Column('info_gid', Integer, ForeignKey('infos.gid'), nullable=False),
    Column('rejecter_id', Integer, ForeignKey('users.user_id'), nullable=False)
)

class Info(Base):
    """Represents a piece of information, a document stored in the system

    Only intrinsic informations that depends only on its content are stored here.
    Access control if done through access_group (who has access) and active_group
    (who trusts the results).
    private_access_group and private_active_group exist to be able to use access groups
    to classify/order the informations, i.e. to still keep track of group ownership even
    for public or no login information.

    A user can also reject a given information adding himself to the rejecters.
    """
    __tablename__='infos'
    #name = models.CharField(max_length=128)
    #description = models.TextField()
    gid = Column(String(57), primary_key=True)
    extraGid = Column(String(130), primary_key=True)
    #internalDirPath = models.CharField(max_length=254)
    #dirPath = models.CharField(max_length=254)
    #immediate_values = models.TextField()
    kind_id = Column(String(57), ForeignKey("info_kinds.gid"))
    kind = relationship("InfoKind",foreign_keys=[kind_id], backref="info_of_kind")
    has_info_of_kind_base = relationship("InfoKind",secondary=info_base_info_kinds,
            backref="in_base_of_infos")
    has_info_of_kind = relationship("InfoKind",secondary=info_info_kinds,
            backref="in_infos")
    #contains = models.ManyToManyField("self", related_name="contained_in", symmetrical=False)
    #related_info = models.ManyToManyField("self")
    #derived_from = models.ManyToManyField("self", related_name="derived_info", symmetrical=False)
    #to_index = models.TextField()
    public = Column(Boolean, default=False)
    access_group_id = Column(Integer, ForeignKey(AccessGroup.access_id))
    access_group = relationship("AccessGroup", foreign_keys=[access_group_id])
    active_group_id = Column(Integer, ForeignKey(AccessGroup.access_id))
    active_group = relationship("AccessGroup", foreign_keys=[active_group_id])
    rejecters = relationship("User", secondary=info_rejecters, backref='rejected_infos')
    #public_override = models.ForeignKey('self', related_name='publicly_overridden', symmetrical=False, null=True, blank=True)

    def kindStr():
        doc = "The kindStr property."
        def fget(self):
            if self.kind:
                return self.kind.name
            else:
                return 'UnknownDocumentType'
        def fset(self, value):
            docs=(KindInfo.objects.filter(name=value,
                kind=basicInfoKinds.gidOf('DocumentType'))
                    .exclude(super_kinds=basicInfoKinds.gidOf('UnknownDocumentType'))).all()
            if len(docs)==0:
                docs=KindInfo.objects.filter(name=value, kind=(basicInfoKinds.gidOf('DocumentType'))).all()
                if not docs:
                    env=InfoKindEnv()
                    kindEl=InfoKindEl(
                        name=value,
                        description='auto generated unknown document type',
                        super_kind='UnknownDocumentType',
                        super_kinds=['UnknownDocumentType'])
                    env.addInfoKindEl(kindEl)
                    env.addDependenciesFrom(basicInfoKinds)
                    KindInfo.ensureInfoKindEnv(env, object_session(self))
                docs=KindInfo.objects.filter(name=value, kind=(basicInfoKinds.gidOf('DocumentType'))).all()
            if docs:
                self.kind=[docs[0]]
        def fdel(self):
            del self._kindStr
        return locals()
    kindStr = property(**kindStr())

    def setKindStr(kindStr, metaTypes):
        """sets kindStr, and if unknown tries to see if metaTypes define it, otherwise
        adds an unknown type for it"""
        docs=(KindInfo.objects.filter(name=kindStr, kind=basicInfoKinds.gidOf('DocumentType'))
        .exclude(super_kinds=basicInfoKinds.gidOf('UnknownDocumentType'))).all()
        if docs:
            if len(docs)>1:
                logging.getLogger("nomad.InfoKind").warn('Multiple KindInfo for Document type %s: %d')
        else:
            metaT=KindInfoEnv(metaTypes)
            if kindStr in metaT:
                # reduce to kindStr + dependencies?
                InfoKind.ensureInfoKindEnv(metaT)
            self.kindStr=kindStr

    def __init__(self, egid=None, kindStr=None,**kwds):
        if not egid is None:
            if egid[-1] != 'I' and egid[-1] != 'T':
                raise Exception('Invalid egid for info '+egid)
            kwds['gid']=egid[:-1]
        super(Info,self).__init__(**kwds)
        if not kindStr is None:
            self.kindStr=kindStr
    def egid():
            doc = "The egid property."
            def fget(self):
                    return self.gid + "I"
            def fset(self, value):
                    if len(value)!= 58 and len(value)!= 34:
                            raise Exception("Invalid length for egid value for info: {0!r}".format(value))
                    if value[-1]!="I":
                        raise Exception("Invalid egid value for info: {0!r}".format(value))
                    self.gid = value[:57]
            return locals()
    egid = property(**egid())

    def update_has_info_of_kind(self):
        newSet=set()
        for k in self.has_info_of_kind_base.all():
            if k.gid in newSet:
                continue
            toDo=[k]
            while (len(toDo)>0):
                kNow=toDo.pop()
                if kNow.gid in newSet:
                    continue
                newSet.add(kNow.gid)
                if kNow.publicOverride is not None and not kNow.publicOverride in newSet:
                    toDo.append(kNow.publicOverride)
                for j in kNow.super_kinds.all():
                    if not j in newSet:
                        toDo.append(j)
        self.has_info_of_kind=newSet

    def visibleUsers():
            doc = "The visibleUsers as space separated list."
            def fget(self):
                return " ".join(self.access_group.users.all().order_by("name").values_list("name", flat=True))
            return locals()
    visibleUsers = property(**visibleUsers())

    def hasAccess(user):
        if user in self.access_group.users.all():
            return True
        return False

    def isActiveForUser(user):
        if user in self.active_group.users.all():
            return True
        return False

    def __unicode__(self):
        return self.name

    def save(self, *argv, **kwds):
        self.update_has_info_of_kind()
        super(InfoKind, self).save(*argv, **kwds)

    def addToGroup(self,accessControlGroup, addToActive=True):
        """Gives access to this information to the given group, and (if addToActive is true)
        accepts also its contents"""
        if not accessControlGroup in self.access_group.control_groups.all():
            groups=list(self.access_group.control_groups.values_list('id').all())
            groups.append(accessControlGroup.id)
            self.access_group=AccessControlGroup.getForGroups(groups)
            if (not accessControlGroup.isPublicOrNoLogin() and
                    not accessControlGroup in self.private_access_group.control_groups.all()):
                groups=list(self.private_access_group.control_groups.values_list('id').all())
                groups.append(accessControlGroup.id)
                self.private_access_group=AccessControlGroup.getForGroups(groups)
        if addToActive and not accessControlGroup in self.active_group.control_groups.all():
            groups=list(self.active_group.control_groups.values_list('id').all())
            groups.append(accessControlGroup.id)
            self.active_group=AccessControlGroup.getForGroups(groups)
            if (not accessControlGroup.isPublicOrNoLogin() and
                    not accessControlGroup in self.private_active_group.control_groups.all()):
                groups=list(self.private_active_group.control_groups.values_list('id').all())
                groups.append(accessControlGroup.id)
                self.private_active_group=AccessControlGroup.getForGroups(groups)

    @staticmethod
    def defaultsForAccessGroupId(controlGroup, baseDict=None):
        """returns a dictionary with the defaults to give access to controlGroup

        if baseDict is given adds the defaults to it.
        Typically used for the defaults argument of get_or_create."""
        if baseDict is None:
            baseDict = {}
        accessGroup=AccessGroup.getForGroups([controlGroup.id])
        if controlGroup.isPublicOrNoLogin():
            emptyGroup=AccessGroup.getForGroups([])
            baseDict['access_group']=accessGroup
            baseDict['private_access_group']=emptyGroup
            baseDict['active_group']=accessGroup
            baseDict['private_active_group']=emptyGroup
        else:
            baseDict['access_group']=accessGroup
            baseDict['private_access_group']=accessGroup
            baseDict['active_group']=accessGroup
            baseDict['private_active_group']=accessGroup
        return baseDict

    def rejectFromGroup(self,accessControlGroup):
        """removes  this information from the given group"""
        if accessControlGroup.isPublicOrNoLogin():
            # FIXME
            assert False
        if accessControlGroup in self.active_group.control_groups.all():
            groups=list(self.active_group.control_groups.values_list('id').all())
            groups.remove(accessControlGroup.id)
            self.active_group=AccessControlGroup.getForGroups(groups)

#class InsertInfo(models.Model):
#  gid = models.CharField(max_length=57, db_index=True)
#  kinds = models.ManyToManyField(InfoKind, related_name='inserters')
#  infos = models.ManyToManyField(Info, related_name='inserters')
#  inserted_by = models.CharField(max_length=254)
#  inserted_on = models.DateField(auto_now_add=True)
#  collected_by = models.CharField(max_length=254)
#  collected_on = models.DateField()
#
#  def __unicode__(self):
#    return u"InsertInfo "+from_info.name+u"->"+to_info.name

class InfoUse(Base):
    __tablename__='info_uses'
    __table_args__=(UniqueConstraint('egid','uri','uri_class',name='u_info_uses1'),)
    id = Column(Integer, primary_key=True)
    egid = Column(String(58),index=True)
    uri = Column(String)
    uri_class = Column(Enum(*[x[1] for x in storage_document.DocClass.choices], name='uri_classes'))
    is_storage = Column(Boolean)

    def __unicode__(self):
        return u"Use "+self.egid+u"->"+self.uri

class Replacers(Base):
    __tablename__="replacers"
    __table_args__=(UniqueConstraint('from_egid','to_egid','defined_in_egid',name='u_replacers1'),)
    id = Column(Integer, primary_key=True)
    from_egid = Column(String(58), index=True)
    to_egid = Column(String(58), index=True)
    defined_in_egid = Column(String(58), index=True)

    def __unicode__(self):
        return u"Replaced "+self.from_gid+u"->"+self.do_gid

#class NormalizedGid(models.Model):
#  from_gid = models.CharField(max_length=58, db_index=True)
#  to_gid = models.CharField(max_length=58, db_index=True)
#  defined_in_gid = models.CharField(max_length=58, db_index=True)
#
#  def __unicode__(self):
#    return u"Replaced "+self.from_gid+u"->"+self.do_gid

#class PurgedInfo(models.Model):
#  access_control_group = models.ForeignKey(AccessControlGroup, related_name='purged_infos',
#    null=True)
#  purgedGid = models.CharField(max_length=58, db_index=True)
#  publicOverride = models.ForeignKey(Info, null=True, blank=True)
#
#  def __unicode__(self):
#    if self.publicOverride is None:
#      return u"Purged "+self.purgedGid
#    else:
#      return u"Purged "+self.purgedGid+u"->"+self.publicOverride.gid

class Connection(Base):
    __tablename__="connections"
    __table_args__=(UniqueConstraint('from_egid','to_egid','defined_in_egid','kind_gid',name='u_connection1'),)
    #access_group = models.ForeignKey(AccessGroup, related_name="connections")
    id = Column(Integer, primary_key=True)
    from_egid = Column(String(58), index=True)
    to_egid = Column(String(58), index=True)
    defined_in_egid = Column(String(58), index=True)
    distance = Column(Float)
    kind_gid = Column(String(57),ForeignKey('info_kinds.gid'))
    kind = relationship('InfoKind', foreign_keys=[kind_gid])

    def __unicode__(self):
        return (u"Connection "+self.kind.name+u" "+self.from_egid+u"->"+self.to_egid+
                u" ("+self.defined_in_egid+u")")
