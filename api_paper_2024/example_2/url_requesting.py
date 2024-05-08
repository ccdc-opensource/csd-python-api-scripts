# -*- coding: utf-8 -*-
#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
"""
A small wrapper to the requests module to implement retries and centralise exceptions
"""
import time
from http import HTTPStatus

import requests
from requests.exceptions import HTTPError


class URLRequestError(RuntimeError):
    def __init__(self, msg, status_code):
        super().__init__(msg)
        self.status_code = status_code


class URLRequest:
    def __init__(self, logger):
        self._logger = logger

    def run(self, url, retries=3):
        retry_codes = [
            HTTPStatus.TOO_MANY_REQUESTS,
            HTTPStatus.INTERNAL_SERVER_ERROR,
            HTTPStatus.BAD_GATEWAY,
            HTTPStatus.SERVICE_UNAVAILABLE,
            HTTPStatus.GATEWAY_TIMEOUT,
        ]
        for n in range(retries):
            try:
                request = requests.get(url)
                request.raise_for_status()
                if self._logger is not None:
                    self._logger.debug(f"{url} - success!")
                break
            except HTTPError as exc:

                code = exc.response.status_code
                if code in retry_codes:
                    if self._logger is not None:
                        self._logger.debug(f"Exception raised {exc} for {url} - will retry")
                    time.sleep(n)
                    continue
                if self._logger is not None:
                    self._logger.warning(f"Unhandleable exception raised {exc} for {url}")
                raise URLRequestError(f"Query failed: {exc} {url}", exc.response.status_code)

        try:
            return request.json()
        except requests.exceptions.JSONDecodeError as exc:
            if self._logger is not None:
                self._logger.warning(f"Decode of request failed {request} {url}")
            raise URLRequestError(f"Decode of request failed {request} {url}", -9999)
