# Overrides #

This directory contains a collection of override files.
These files describe the new version of a NomadMetaInfo, by listing old gid and new gid, and can thus introduce versioning for NomadMetaInfo.
The name or other keys can be given, but are only informative and can be omitted.

## File Naming convention ##

Normally overrides are given between two tagged versions or between the last checked in version,
and the current state. so override files are by default given as

    <oldVersion>_<newVersion>_YYYY-MM-DD.nomadmetainfo_overrides.json

where *oldVersion* can be the first 10 characters of the git sha, a tag name, or even empty;
just like *newVersion*. *YYYY-MM-DD* is the current date, and if required an "\_n" with a suitable
number *n* that does not clash with existing files can be used.

The extension .nomadmetainfo_overrides.json is mandatory.

## Automatic Generation ##

Normally you can generate these files automatically with scripts/nomadscripts/calculate_meta_info_overrides.py
The script works if the name of the NomadMetaInfo is the same but have differen gid.

## Complex Cases ##

In cases in which you have renamed a NomadMetaInfo or there is a NomadMetaInfo outside the
standard that you want to replace with the standard one you have to create (or complete)
the override file by hand.
In these cases the output with --verbose can be useful.
It is also possible to use scripts/nomadscripts/normalize_local_kinds.py --add-gid to update
each NomadMetaInfo with its gid, which then you can use to create manual override files.
Please do not check in the repository the generated .nomadmetainfo.json files with gid.
