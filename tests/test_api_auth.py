import sys
import os
import unittest
import importlib

# We need to test the Flask app in isolation without it loading a real GFA file.
# Patch the data path before importing server.py so no FileNotFoundError occurs.

# Stub out the GFA load by monkeypatching BipartiteGraph before import
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'viz'))

# --------------------------------------------------------------------------
# We need to prevent server.py from loading a real GFA at import time.
# We do this by providing a stub BipartiteGraph that does nothing on parse_gfa.
# --------------------------------------------------------------------------
import bipartite as _bipartite_real
_real_bipartite_graph = _bipartite_real.BipartiteGraph


class _StubGraph:
    states = {}

    def parse_gfa(self, path):
        pass

    def compute_addresses(self):
        pass


_bipartite_real.BipartiteGraph = lambda engine, store: _StubGraph()

# Now safe to import server
import importlib.util
_spec = importlib.util.spec_from_file_location(
    "server",
    os.path.join(os.path.dirname(__file__), '..', 'viz', 'server.py'),
)
_server_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_server_mod)
app = _server_mod.app
_bipartite_real.BipartiteGraph = _real_bipartite_graph


# ======================================================================
# Helpers
# ======================================================================

def _client(token: str = ""):
    """Return a Flask test client with FRX_API_TOKEN set to token."""
    orig = _server_mod._API_TOKEN
    _server_mod._API_TOKEN = token
    client = app.test_client()
    return client, orig


# ======================================================================
# /health - always public
# ======================================================================

class TestHealthEndpoint(unittest.TestCase):
    """Health endpoint must be public regardless of auth config."""

    def setUp(self):
        self.client = app.test_client()

    def test_health_200_no_auth(self):
        _server_mod._API_TOKEN = "secret"
        resp = self.client.get('/health')
        self.assertEqual(resp.status_code, 200)

    def test_health_200_dev_mode(self):
        _server_mod._API_TOKEN = ""
        resp = self.client.get('/health')
        self.assertEqual(resp.status_code, 200)

    def test_health_returns_json(self):
        _server_mod._API_TOKEN = ""
        resp = self.client.get('/health')
        data = resp.get_json()
        self.assertIn('status', data)
        self.assertEqual(data['status'], 'ok')

    def test_health_auth_field_true_when_token_set(self):
        _server_mod._API_TOKEN = "secret"
        resp = self.client.get('/health')
        data = resp.get_json()
        self.assertTrue(data['auth'])

    def test_health_auth_field_false_in_dev_mode(self):
        _server_mod._API_TOKEN = ""
        resp = self.client.get('/health')
        data = resp.get_json()
        self.assertFalse(data['auth'])


# ======================================================================
# /graph - auth enforcement
# ======================================================================

class TestGraphAuth(unittest.TestCase):
    """Graph endpoint must enforce Bearer token auth when configured."""

    def setUp(self):
        self.client = app.test_client()

    def tearDown(self):
        _server_mod._API_TOKEN = ""

    def test_graph_200_dev_mode(self):
        _server_mod._API_TOKEN = ""
        resp = self.client.get('/graph')
        self.assertEqual(resp.status_code, 200)

    def test_graph_401_no_header_when_auth_enabled(self):
        _server_mod._API_TOKEN = "correct_token"
        resp = self.client.get('/graph')
        self.assertEqual(resp.status_code, 401)

    def test_graph_401_wrong_token(self):
        _server_mod._API_TOKEN = "correct_token"
        resp = self.client.get('/graph', headers={'Authorization': 'Bearer wrong_token'})
        self.assertEqual(resp.status_code, 401)

    def test_graph_200_correct_token(self):
        _server_mod._API_TOKEN = "correct_token"
        resp = self.client.get('/graph', headers={'Authorization': 'Bearer correct_token'})
        self.assertEqual(resp.status_code, 200)

    def test_graph_401_malformed_header(self):
        _server_mod._API_TOKEN = "correct_token"
        # Missing "Bearer " prefix
        resp = self.client.get('/graph', headers={'Authorization': 'correct_token'})
        self.assertEqual(resp.status_code, 401)

    def test_graph_401_empty_token(self):
        _server_mod._API_TOKEN = "correct_token"
        resp = self.client.get('/graph', headers={'Authorization': 'Bearer '})
        self.assertEqual(resp.status_code, 401)

    def test_graph_401_returns_json_error(self):
        _server_mod._API_TOKEN = "correct_token"
        resp = self.client.get('/graph')
        data = resp.get_json()
        self.assertIn('error', data)
        self.assertEqual(data['error'], 'Unauthorized')


# ======================================================================
# /query - auth enforcement
# ======================================================================

class TestQueryAuth(unittest.TestCase):
    """Query endpoint must enforce the same auth as /graph."""

    def setUp(self):
        self.client = app.test_client()

    def tearDown(self):
        _server_mod._API_TOKEN = ""

    def test_query_401_no_auth(self):
        _server_mod._API_TOKEN = "secret"
        resp = self.client.get('/query?path=Root/1')
        self.assertEqual(resp.status_code, 401)

    def test_query_allowed_with_correct_token(self):
        _server_mod._API_TOKEN = "secret"
        resp = self.client.get('/query?path=Root/1',
                               headers={'Authorization': 'Bearer secret'})
        # May be 200 or 400 (no results) but NOT 401
        self.assertNotEqual(resp.status_code, 401)

    def test_query_400_no_params_but_authenticated(self):
        _server_mod._API_TOKEN = ""
        resp = self.client.get('/query')
        self.assertEqual(resp.status_code, 400)


# ======================================================================
# Rate limiting
# ======================================================================

class TestRateLimiting(unittest.TestCase):
    """Rate limiter must return 429 after the limit is breached."""

    def setUp(self):
        self.client = app.test_client()
        _server_mod._API_TOKEN = ""
        # Clear rate registry to start fresh
        _server_mod._rate_registry.clear()

    def tearDown(self):
        _server_mod._rate_registry.clear()
        _server_mod._API_TOKEN = ""

    def test_first_request_not_rate_limited(self):
        orig = _server_mod._RATE_LIMIT
        _server_mod._RATE_LIMIT = 5
        try:
            resp = self.client.get('/health')
            self.assertNotEqual(resp.status_code, 429)
        finally:
            _server_mod._RATE_LIMIT = orig

    def test_rate_limit_429_after_limit(self):
        orig = _server_mod._RATE_LIMIT
        _server_mod._RATE_LIMIT = 3
        _server_mod._rate_registry.clear()
        try:
            # /health is NOT rate-limited; use /graph instead
            for _ in range(3):
                self.client.get('/graph')
            # 4th should be rate-limited
            resp = self.client.get('/graph')
            self.assertEqual(resp.status_code, 429)
        finally:
            _server_mod._RATE_LIMIT = orig
            _server_mod._rate_registry.clear()

    def test_rate_limit_response_has_retry_after(self):
        orig = _server_mod._RATE_LIMIT
        _server_mod._RATE_LIMIT = 1
        _server_mod._rate_registry.clear()
        try:
            self.client.get('/graph')
            resp = self.client.get('/graph')
            if resp.status_code == 429:
                self.assertIn('Retry-After', resp.headers)
        finally:
            _server_mod._RATE_LIMIT = orig
            _server_mod._rate_registry.clear()

    def test_rate_limit_returns_json(self):
        orig = _server_mod._RATE_LIMIT
        _server_mod._RATE_LIMIT = 1
        _server_mod._rate_registry.clear()
        try:
            self.client.get('/graph')
            resp = self.client.get('/graph')
            if resp.status_code == 429:
                data = resp.get_json()
                self.assertIn('error', data)
                self.assertEqual(data['error'], 'TooManyRequests')
        finally:
            _server_mod._RATE_LIMIT = orig
            _server_mod._rate_registry.clear()

    def test_health_not_rate_limited(self):
        orig = _server_mod._RATE_LIMIT
        _server_mod._RATE_LIMIT = 2
        _server_mod._rate_registry.clear()
        try:
            # Hit graph limit
            for _ in range(3):
                self.client.get('/graph')
            # /health should still be 200
            resp = self.client.get('/health')
            self.assertEqual(resp.status_code, 200)
        finally:
            _server_mod._RATE_LIMIT = orig
            _server_mod._rate_registry.clear()


if __name__ == '__main__':
    unittest.main(verbosity=2)
