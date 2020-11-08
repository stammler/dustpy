import json
import requests
import socket
from simframe.utils import colorize


def print_version_warning(timeout=0.5):
    """Functions prints a warning to screen if a newer version of ``DustPy`` is available.

    Parameters
    ----------
    timeout : float, optional, default : 0.5
        Timeout in seconds after which test will be aborted

    Notes
    -----
    If no internet connection is available the test will be aborted after ``timeout``.
    This also happens if an internet connection is available but the ping is greater than ``timeout``."""

    # Only check for version if internet connection is available
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect(
            ("8.8.8.8", 53))

        from dustpy.utils import __version__ as ver_cur

        data = requests.get("https://pypi.org/pypi/dustpy/json").text
        ver_lat = json.loads(data)["info"]["version"]

        if ver_cur != ver_lat:
            msg = "\nA newer version of DustPy is available.\nThis version:   {}\nLatest version: {}\n\nUpgrade with\npip install dustpy --upgrade\n".format(
                ver_cur, ver_lat)
            msg = colorize(msg, color="purple")
            print(msg)
    except:
        pass
