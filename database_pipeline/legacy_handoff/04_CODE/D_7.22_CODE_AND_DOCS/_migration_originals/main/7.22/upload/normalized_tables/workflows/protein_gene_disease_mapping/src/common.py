from __future__ import annotations

import csv
import hashlib
import json
import logging
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import yaml


def find_repo_root(start: Path | None = None) -> Path:
    current = (start or Path(__file__)).resolve()
    if current.is_file():
        current = current.parent
    for candidate in (current, *current.parents):
        if (candidate / ".git").exists():
            return candidate
    raise RuntimeError("Repository root (.git) not found")


def load_config(path: str | Path) -> tuple[Path, dict[str, Any]]:
    root = find_repo_root(Path(path))
    config_path = Path(path)
    if not config_path.is_absolute():
        config_path = root / config_path
    data = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    return root, data


def repo_path(root: Path, value: str | Path) -> Path:
    p = Path(value)
    return p if p.is_absolute() else root / p


def ensure_dirs(config: dict[str, Any], root: Path) -> dict[str, Path]:
    workflow = repo_path(root, config["outputs"]["workflow_dir"])
    release = repo_path(root, config["outputs"]["release_dir"])
    paths = {
        "workflow": workflow,
        "release": release,
        "cache": workflow / "cache",
        "interim": workflow / "interim",
        "logs": workflow / "logs",
        "reports": workflow / "reports",
    }
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    return paths


def configure_logging(log_path: Path) -> logging.Logger:
    logger = logging.getLogger("protein_gene_disease_mapping")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s")
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    logger.addHandler(stream)
    return logger


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def retrieval_date() -> str:
    return datetime.now(timezone.utc).date().isoformat()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        materialized = list(rows)
        fields = list(materialized[0]) if materialized else []
        row_iterable: Iterable[dict[str, Any]] = materialized
    else:
        row_iterable = rows
    with path.open("w", encoding="utf-8", newline="") as handle:
        if not fields:
            return
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in row_iterable:
            writer.writerow({key: clean_cell(row.get(key, "")) for key in fields})


def clean_cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return str(value).replace("\x00", "").replace("\r\n", " ").replace("\r", " ").replace("\n", " ")


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_json(path: Path, default: Any = None) -> Any:
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_accession(accession: str) -> str:
    return re.sub(r"-\d+$", "", accession.strip())


def split_multi(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v).strip() for v in value if str(v).strip()]
    return [x.strip() for x in re.split(r"[;|,]", str(value)) if x.strip()]


def stable_join(values: Iterable[Any], sep: str = ";") -> str:
    unique = sorted({str(v).strip() for v in values if str(v).strip()})
    return sep.join(unique)


def strict_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (value or "").lower()).strip()


def request_with_retry(session, method: str, url: str, *, max_retries: int, timeout: int,
                       backoff: float, logger: logging.Logger, **kwargs):
    last_error: Exception | None = None
    for attempt in range(max_retries):
        try:
            response = session.request(method, url, timeout=timeout, **kwargs)
            if response.status_code in {429, 500, 502, 503, 504}:
                raise RuntimeError(f"HTTP {response.status_code}: {response.text[:300]}")
            if 400 <= response.status_code < 500:
                response.raise_for_status()
            response.raise_for_status()
            return response
        except Exception as exc:  # requests exposes several transport exceptions
            last_error = exc
            status_code = getattr(getattr(exc, "response", None), "status_code", None)
            if status_code is not None and 400 <= status_code < 500 and status_code not in {400, 408, 409, 429}:
                break
            if attempt + 1 == max_retries:
                break
            delay = backoff * (2 ** attempt)
            logger.warning("request retry %s/%s for %s after %s", attempt + 1, max_retries, url, exc)
            time.sleep(delay)
    raise RuntimeError(f"Request failed after {max_retries} attempts: {url}: {last_error}")


def nonempty(value: Any) -> bool:
    return bool(str(value or "").strip())
