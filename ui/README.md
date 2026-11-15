# VariantNavigator local UI

Run from the repository root:

```powershell
python ui/server.py --port 8000
```

Open http://127.0.0.1:8000. Opening index.html directly cannot call the API.

## Real pipeline behavior

- Validation reads `output/model/validation_predictions.json` on each request. Metrics, scatterplot, table and CSV export use this artifact. Only its saved head has per-row estimates; other heads are not fabricated.
- Submissions accept one standard amino-acid substitution at ABL1 positions 242?512. The pipeline verifies reference identity.
- Known variants return real measured net growth. An independently saved held-out estimate is shown alongside it when available.
- Unmeasured variants invoke `pipeline.predict.predict_variant`, which fits and caches the existing handcrafted-feature RandomForest. It requires NumPy/scikit-learn and may access public structure/sequence resources through pipeline feature extraction. The first fit can take time.
- Runtime RandomForest estimates do not inherit saved ESM validation metrics. This server does not launch Modal/GPU jobs or use the persisted ESM head for novel inference.
- Run overview shows the validation artifact SHA-256 and its provenance limitations. Replicate intervals are not confidence intervals.

For runtime inference, install into the Python environment used to run the server:

```powershell
python -m pip install numpy scikit-learn
```

The UI uses native fetch and a Python standard-library loopback server. It serves only the ui directory; credentials and model files are not static assets. Requests are serialized during prediction to protect cached model fitting. There is no public hosting or authentication; keep the server bound to loopback.

## Verification

```powershell
python -m unittest discover -s tests -p test_ui_server.py
node --check ui/app.js
```

The poster theme and helix decorations are retained. No synthetic model responses remain.
