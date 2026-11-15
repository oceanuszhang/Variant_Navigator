"""API contract tests without network, GPU, or model fitting."""
import json
from threading import Thread
import unittest
from unittest.mock import patch
from urllib.error import HTTPError
from urllib.request import Request, urlopen
from http.server import ThreadingHTTPServer

from ui.server import Handler, submit


class UIServerTests(unittest.TestCase):
    def setUp(self):
        self.server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        self.thread = Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.url = f'http://127.0.0.1:{self.server.server_port}'

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()

    def test_validation_and_static(self):
        with patch('ui.server.validation', return_value={'rows': [], 'head': 'ridge'}):
            with urlopen(self.url+'/api/validation') as response:
                self.assertEqual(json.load(response)['head'], 'ridge')
        with urlopen(self.url+'/') as response:
            self.assertIn(b'VariantNavigator', response.read())

    def test_submission_routes_to_pipeline(self):
        with patch('ui.server.submit', return_value={'result': {'kind': 'inference'}}) as inference:
            request = Request(self.url+'/api/predict', data=b'{"variant":" f486s ","sequence":"ABC"}', headers={'Content-Type':'application/json'})
            with urlopen(request) as response:
                self.assertEqual(json.load(response)['result']['kind'], 'inference')
            inference.assert_called_once_with('F486S', 'ABC')

    def test_invalid_json_and_cross_origin(self):
        for body, headers, expected in [(b'[]', {}, 400), (b'{}', {'Origin':'https://example.com'}, 403)]:
            with self.assertRaises(HTTPError) as error:
                urlopen(Request(self.url+'/api/predict', data=body, headers=headers))
            self.assertEqual(error.exception.code, expected)
            error.exception.close()

    def test_invalid_variant_rejected(self):
        with self.assertRaises(ValueError):
            submit('../secrets')
        with self.assertRaisesRegex(ValueError, 'is T, not A'):
            submit('A315I')
        with self.assertRaisesRegex(ValueError, 'study reference'):
            submit('T315I', 'ACDE')
        with self.assertRaisesRegex(ValueError, 'standard residue'):
            submit('T315I', 'ABC123')

    def test_measurement_separate_from_held_out(self):
        with patch('pipeline.predict.predict_variant', return_value={'kind':'measurement', 'measured_netgrowth':0.04}), patch('pipeline.predict.held_out', return_value={'predicted':0.01}):
            result = submit('T315I')
            self.assertEqual(result['result']['kind'], 'measurement')
            self.assertEqual(result['held_out']['predicted'], 0.01)
            self.assertEqual(result['sequence']['ref'], 'T')
            from pipeline import predict
            fasta = '>ABL1\n' + predict.REFERENCE.read_text().strip().lower()
            self.assertEqual(submit('T315I', fasta)['result']['kind'], 'measurement')
