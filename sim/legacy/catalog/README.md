# Field registry

Each site owns its survey references. Each acquisition owns an explicit check
design and output report path. Coordinates and dimensions carry their datums;
missing target masks and field detections remain JSON `null`.

`surveys.json` geometry is a **historical reference**, not data inferred from SAR.
Models are simplified prisms in local metres. Every model requires a primary
source, a vertical datum, disclosed assumptions and acquisition-state notes.
Additional candidates can be listed without constructing a model.

Do not put raw acquisitions, PDF archives, local user paths or caches here.
The Giza landmark record contains selected published coordinates and a tiny
geoid-grid extract, with source hashes and URLs. Keep archive downloads and
large working data outside the repository.

See [the field workflow](../fieldwork/README.md) for reproducible commands,
schema conventions, current limitations and how to add another site or scan.
