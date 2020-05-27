#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
@package    nomadcore.nomad_query
@copyright  Copyright 2018+ Fritz Haber Institute of the Max Planck Society,
            Benjamin Regler - Apache 2.0 License
@license    http://www.apache.org/licenses/LICENSE-2.0
@author     Benjamin Regler
@version    1.0.0

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import os
import re
import sys
import json
import time
import errno
import random

if sys.version_info.major > 2:
    # For Python 3.0 and later
    from urllib.parse import quote, unquote_plus
    from urllib.request import urlopen, Request
else:
    # Fall back to Python 2's urllib2
    from urllib import quote, unquote_plus
    from urllib2 import urlopen, Request


class NomadQueryResult(object):
    """Nomad Query Result class

    @author Benjamin Regler
    """

    def __init__(self, query, response, version=1.0):
        """Constructor.
        
        Arguments:
            query {dict}    -- Query information, i.e., query filter, context,
                               group_by, endpoint, and URL
            response {dict} -- Response of the Nomad Query API
        
        Keyword Arguments:
            version {number} -- Version of the Nomad Query data file
                                (default: {1.0})
        """
        self._uri = []
        self._query = query or {}
        self._timestamp = int(time.time())
        self._response = response.get('result', {})
        self._version = version

        # Construct download path
        path = response.get('path', '')
        self._download_url = self._query.get('endpoint', '') + 'download/' + \
            path.split('_')[-1] + '?file=' + quote(path.encode('utf-8')) + '.json'

        # Get Nomad URIs
        response = NomadQuery().request(self._download_url)
        if response['status'] == 'success':
            regex = re.compile(r'(?<=/[a-zA-Z0-9\-_]{3}/)[^\.]+')
            paths = response['data'].get('result', [])

            for path in paths:
                match = regex.search(path)
                if match:
                    # Substitute prefixes
                    groups = match.group(0).split('/')
                    groups[0] = 'N' + groups[0][1:]         # Normalized

                    if len(groups) == 2:
                        groups[1] = 'C' + groups[1][1:]     # Computed

                self._uri.append('nmd://' + '/'.join(groups))

    def version(self):
        """Get the version of the Nomad Query data file.
        
        Returns:
            float -- Version of the Nomad Query data file
        """
        return self._version

    def timestamp(self):
        """Get the timestamp of the query.

        Returns:
            tuple -- The timestamp of the query
        """
        return self._timestamp

    def download_url(self):
        """Get the download URL of the query.

        Returns:
            str -- The download URL of the query
        """
        return self._download_url

    def query(self):
        """Get query information

        Returns:
            dict -- Query information
        """
        return self._query

    def response(self):
        """Get the query response.

        Returns:
            dict -- The query response
        """
        return self._response

    def uri(self):
        """Get Nomad URIs.

        Returns:
            list -- List of Nomad URIs
        """
        return self._uri


class NomadQuery(object):
    """Nomad Query class for accessing stored queries

    @author Benjamin Regler
    """

    # Version of the Nomad Query API
    __version__ = 1.0

    # Nomad API endpoint
    endpoint = os.environ.get('NOMAD_BASE_URI','https://analytics-toolkit.nomad-coe.eu') + '/api/'

    # Private user path
    user_path = '/data/private'

    def __init__(self, username='', endpoint=''):
        """Constructor.

        Keyword Arguments:
            username {str} -- Current username. Leave empty to auto-detect
                              username (default: {''})
            endpoint {str} -- Endpoint of the Nomad API (default:
                              ${NOMAD_BASE_URI}/api if set, otherwise
                              {'https://analytics-toolkit.nomad-coe.eu/api/'})
        """
        self._username = ''
        self._base_path = ''

        # Guess username (more like a HACK)
        if not username:
            if os.path.exists(self.user_path):
                paths = os.listdir(self.user_path)
                if len(paths) == 1 and paths[0].lower() != 'nomad':
                    username = paths[0]

        # Set username and overwrite endpoint, if required
        self.username(username)
        if endpoint:
            self.endpoint = str(endpoint)

    def username(self, username=''):
        """Get or set the username.

        Keyword Arguments:
            username {str} -- Current username (default: {''})

        Returns:
            str -- The current username
        """
        if username:
            self._username = str(username)
            self._base_path = os.path.join(self.user_path, self._username,
                                           'nomad-query')
        return self._username

    def request(self, url, timeout=10):
        """Request a URL

        Arguments:
            url {str} -- The URL of a web address

        Keyword Arguments:
            timeout {number} -- Timeout of the request in seconds (default: {10})

        Returns:
            dict -- A dictionary with success status, response data, or 
                    error message
        """
        # Default request response
        result = {
            'url': url,
            'status': 'error',
            'message': 'Unknown error. Please inform the Nomad team to '
                       'solve this problem.'
        }

        try:
            # Get URL
            response = urlopen(Request(url), timeout=timeout)
            if response.code != 200:
                raise RuntimeError(result['message'])

            # Read response
            data = json.loads(response.read().decode('utf-8'), 'utf-8')

            # Populate result
            result.pop('message')
            result.update({
                'status': 'success',
                'data': data
            })
        except Exception as exc:
            exc = sys.exc_info()[1]
            response = result.copy()

            # Get error message
            message = exc
            if sys.version_info <= (2, 5) and hasattr(exc, 'message'):
                message = exc.message
            elif hasattr(exc, 'reason'):
                message = exc.reason
            response['message'] = str(message)

            # Fix error message
            if response['message'].endswith('timed out'):
                response['message'] = 'Connection timed out. The Nomad ' + \
                    'Analytics API Service is currently unavailable.'

        # Return result
        return result

    def resolve(self, nmd, path='', recursive=True, timeout=10):
        """Resolve a Nomad URI.

        A Nomad URI is a base32 encoded string usually 28-characters wide long.
        If separated with a slash, the first URI represents the repository and
        the second URI represents the calculation, e.g.,

            nmd://N2gHWfGLtvdPhqTq8jL6wq1GeLV7_/C093CpZoWf8U32dmRFkp17gy4BCvT .

        Arguments:
            nmd {str} -- Nomad URI with prefix `nmd://`

        Keyword Arguments:
            path {str}       -- Path for resolving Nomad URI, i.e.,
                                /section_run/0c/section_single_configuration_calculation/0c
                                (default: {''})
            recursive {bool} -- Recursively resolve Nomad URI. This may take
                                longer and uses much more memory
                                (default: {True})
            timeout {number} -- Timeout of the request in seconds (default: {10})

        Returns:
            dict -- A python dictionary with data for the given Nomad URI

        Raises:
            ValueError -- Invalid scheme. Nomad URI must start with "nmd://".
            RuntimeError -- Resolving Nomad URI "..." failed.
        """
        if not nmd.startswith('nmd://'):
            raise ValueError('Invalid scheme. Nomad URI must start with "nmd://".')

        # Construct URL
        url = self.endpoint + 'resolve/' + nmd[6:].strip('/')
        if path:
            url += '/' + path.strip('/')
        if recursive:
            url += '?format=recursiveJson'

        response = self.request(url, timeout=timeout)
        if response['status'] != 'success':
            raise RuntimeError('Resolving Nomad URI "%s" failed.' % nmd)

        return response['data']

    def list(self):
        """List all stored queries.

        Returns:
            list -- A sorted list of stored queries. Index 0 always refers to
                    the latest query. Index 1 to the second latest and so on
                    and so forth
        """
        queries = []
        base_path = os.path.join(self._base_path, 'data')
        if not os.path.isdir(base_path):
            return queries

         # Get all stored queries
        for filename in os.listdir(base_path):
            path = os.path.join(base_path, filename)
            if os.path.isfile(path):
                modified = os.path.getmtime(path)
                name = os.path.splitext(filename)[0].split(' ')[-1]

                # Store some useful information about the query
                queries.append({
                    'name': name,
                    'path': path,
                    'filename': filename,
                    'timestamp': int(modified)
                })

        # Sort queries based on modification time
        queries.sort(key=lambda x: -x['timestamp'])
        return queries

    def query(self, query, group_by='', context='', timeout=10):
        """Query the Nomad Database.

        Arguments:
            query {str} -- The query string (see Nomad API reference)

        Keyword Arguments:
            group_by {str}   -- Group-by field. (default: {''})
            context {str}    -- Query context. Leave empty to use
                                `single_configuration_calculation` (default: {''})
            timeout {number} -- Timeout of the request in seconds (default: {10})

        Returns:
            NomadQueryResult -- The Nomad query result

        Raises:
            RuntimeError -- Connection timed out. The Nomad Analytics API
                            Service is currently unavailable.
            RuntimeError -- Unknown error. Please inform the Nomad team to
                            solve this problem.
        """
        # Set default context
        if not context:
            context = 'single_configuration_calculation'

        # Construct URL
        url = self.endpoint + ('queryGroup/' if group_by else 'query/') + context

        # Add query
        url += '?filter=' + quote(query.strip())
        if group_by:
            url += quote(' GROUPBY ' + group_by.strip().lower())

        # Read URL
        response = self.request(url, timeout=timeout)
        if response['status'] != 'success':
            raise RuntimeError(response['message'])

        # Check connection timeout
        response = response['data']
        if 'timed_out' in response['result'] and response['result']['timed_out']:
            response['message'] = 'Connection timed out.'

        # Check for additional error messages
        if 'message' in response or 'msg' in response:
            raise RuntimeError(response.get('message', response['msg']))

        # Construct Nomad Query response
        query = {
            'context': context,
            'endpoint': self.endpoint,
            'filter': query.strip(),
            'group_by': group_by.strip().lower(),
            'url': url
        }
        return NomadQueryResult(query, response, self.__version__)

    def fetch(self, name_or_index='', resolve=False, **params):
        """Fetch stored query.

        Keyword Arguments:
            name_or_index {int|str} -- The name or the index of the query.
                                       Leave empty to use the latest one
                                       (default: {''})
            resolve {bool}          -- Automatically resolve Nomad URIs
                                       (default: {False})
            params {dict}           -- Parameter passed to resolve method. May
                                       include `size` and `seed` for fetching
                                       only a limited number of URIs chosen by
                                       a given seed (default: {{}})

        Returns:
            dict|bool -- A dictionary with query information, Nomad URIs, and
                         optionally resolved Nomad URIs or False if not found

        Raises:
            KeyError -- Query with name "..." does not exists.
        """
        # Get all queries
        queries = self.list()
        if not queries:
            return False

        # Empty name refers to latest query
        index = name_or_index
        if not name_or_index:
            index = 0

        # Extract path
        path = ''
        if isinstance(index, int):
            # Resolve query index
            path = queries[index].get('path')
        else:
            # Load query with specified name
            for i, query in enumerate(queries):
                if query['name'] == index:
                    path = query['path']
                    index = i
                    break

        # Check if path really exists
        if not os.path.exists(path):
            raise KeyError('Query with name "%s" does not exists.' % name_or_index)

        # Load query - JSON
        with open(path, 'r') as file:
            query = json.load(file, 'utf-8')

        # Add useful variables to query
        query.update({
            'data': {},
            'path': path,
            'name': queries[index].get('name')
        })

        # Resolve Nomad URIs?
        if resolve:
            query['data'] = self._resolve(query['uri'], **params)
        return query

    def save(self, name, query, resolve=False, **params):
        """Save query.

        Arguments:
            name {str}     -- Query name
            query {object} -- Instance of `NomadQueryResult`

        Keyword Arguments:
            resolve {bool}  -- Automatically resolve Nomad URIs
                               (default: {False})
            params {dict}   -- Parameter passed to resolve method. May include
                               `size` and `seed` for fetching only a limited
                               number of URIs chosen by a given seed (default: {{}})

        Returns:
            dict -- A dictionary with query information, Nomad URIs, and
                    optionally resolved data

        Raises:
            RuntimeError -- Query name must be an alphanumeric string and may
                            optionally contain underscores and hyphens.
            RuntimeError -- Query parameter must be a `NomadQueryResult` instance.
        """
        # Check query name
        regex = re.compile(r'[^0-9a-z_\-]+', re.I)
        if not name or regex.search(name):
            raise RuntimeError('Query name must be an alphanumeric string and '
                               'may optionally contain underscores and hyphens.')

        # Check query result
        if not isinstance(query, NomadQueryResult):
            raise RuntimeError('Query parameter must be a `NomadQueryResult` '
                               'instance.')

        # Collect information
        timestamp = query.timestamp()
        now = time.localtime(timestamp)
        filename = time.strftime('%Y-%m-%d %H%M%S ', now) + name + '.json'

        data = {
            'filename': filename,
            'version': self.__version__,
            'query': query.query(),
            'timestamp': timestamp,
            'uri': query.uri()
        }

        # Save file
        path = os.path.join(self._base_path, 'data', filename)

        # Write query result
        self._makedir(os.path.dirname(path))
        with open(path, 'w') as file:
            json.dump(data, file, sort_keys=True, indent=4)

        # Add useful variables to query
        data.update({
            'data': {},
            'path': path,
            'name': name
        })

        # Adjust modification time
        os.utime(path, (timestamp, timestamp))

        # Write log file
        path = os.path.join(self._base_path, 'nomad-query.log')
        url = query.query().get('url', '')

        # Create file and check log size
        if not os.path.exists(path):
            self._makedir(os.path.dirname(self.endpoint.rstrip('/')))

            with open(path, 'w') as file:
                file.write('# This file is auto-generated by the Nomad ' +
                           'Query GUI\n# ' + self.endpoint + '\n\n')

        # Append file
        if url:
            with open(path, 'a') as file:
                file.write(time.strftime('[%Y-%m-%d %H:%M:%S] ', now) +
                           unquote_plus(url) + '\n')

        # Resolve Nomad URIs?
        if resolve:
            data['data'] = self._resolve(data['uri'], **params)
        return data

    def _resolve(self, paths, size=None, seed=None, **params):
        """[Protected] Resolve Nomad URIs.

        Arguments:
            paths {list}  -- List of Nomad URIs

        Keyword Arguments:
            size {number} -- Total number of URIs to resolve (default: {None})
            seed {number} -- Seed to initialize the internal state of the
                             random number generator for selecting n-th Nomad
                             URIs (default: {None})
            params {dict} -- Parameter passed to resolve method (default: {{}})

        Returns:
            dict -- A dictionary with resolved Nomad URIs
        """
        size = params.pop('size', size)
        seed = params.pop('seed', seed)

        # Truncate number of URIs to specified size
        if size and size > 1:
            state = random.getstate()

            # Set seed of random number generator
            if seed is not None:
                random.seed(seed)

            # Get limited set of Nomad URIs
            paths = random.sample(paths, size)
            random.setstate(state)

        return {path: self.resolve(path, **params) for path in paths}

    def _makedir(self, path):
        """[Protected] Recursively create directory.

        This function recursively creates the directory path without raising an
        exception if it does not exists.

        Arguments:
            path {str} -- The path of the directory

        Returns:
            bool -- Returns True on success, False otherwise
        """
        try:
            os.makedirs(path)
        except OSError as exc:  # Python > 2.5
            if exc.errno != errno.EEXIST or not os.path.isdir(path):
                return False
        return True
