"""Normalize locally acquired caption evidence; never download or publish video.

Raw json3 caption files are fetched beforehand with yt-dlp into the ignored
work/acquisition/raw folder (see WORKFLOW.md). This script normalizes them to
timed text and writes the public source manifest with digests.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import gzip
import json
import re

ROOT = Path(__file__).resolve().parents[1]
PRIVATE = ROOT / "work" / "acquisition"
PLAYLIST = ("https://www.youtube.com/playlist?"
            "list=PLLy_2iUCG87ChONRU6EeO6A2rIzM3HUvG")


# ponytail: one hard-coded supplement; read the NPTEL outline if more appear.
SUPPLEMENT = {"id": "NYo2wbePAnI", "duration": 2091,
              "title": "Lecture 44| Inner Product Spaces  - I |Prof Premananda "
                       "Bera | IIT Roorkee | NPTEL"}


def lecture_number(title):
    match = re.search(r"Lecture\s*0*(\d+)", title)
    return int(match[1]) if match else None


def normalize_json3(raw):
    """Keep timed events in order; join their text fragments."""
    result = []
    for event in raw.get("events", []):
        text = "".join(s.get("utf8", "") for s in event.get("segs", []))
        text = re.sub(r"\s+", " ", text).strip()
        if text:
            result.append({"start": event.get("tStartMs", 0) / 1000,
                           "duration": event.get("dDurationMs", 0) / 1000,
                           "text": text})
    if not result:
        raise ValueError("Empty timed caption JSON")
    return result


def entry_record(position, entry, old):
    vid = entry["id"]
    title = entry["title"]
    meta = dict(position=position, playlist_position=position, video_id=vid,
                title=title, duration=entry.get("duration"),
                duration_seconds=entry.get("duration"),
                url=f"https://www.youtube.com/watch?v={vid}",
                lecture_number=lecture_number(title),
                role="introduction" if position == 1 else "lecture")
    raws = sorted((PRIVATE / "raw").glob(f"{position:02d}-{vid}.*json3"))
    if not raws:
        meta.update(status="unavailable")
        return meta
    payload = raws[0].read_bytes()
    track = raws[0].name.split(".")[-2]
    segments = normalize_json3(json.loads(payload.decode("utf-8")))
    archive = PRIVATE / f"{position:02d}-{vid}.json.gz"
    archive.write_bytes(gzip.compress(payload, mtime=0))
    text = "\n".join(f'[{s["start"]:.2f}] {s["text"]}'
                     for s in segments) + "\n"
    normalized = PRIVATE / f"{position:02d}-{vid}.txt"
    normalized.write_bytes(text.encode("utf-8"))
    meta.update(status="acquired", language="en",
                is_generated=not track.startswith("en-"),
                track_name=("English - NPTEL Verified (manual YouTube track "
                            f"{track})" if track.startswith("en-") else
                            "YouTube automatic English captions"),
                acquisition_method="yt-dlp json3 caption download, "
                                   "no video",
                segment_count=len(segments),
                first_start=segments[0]["start"],
                last_end=max(s["start"] + s["duration"] for s in segments),
                raw_sha256=sha256(payload).hexdigest(),
                normalized_text_sha256=sha256(
                    normalized.read_bytes()).hexdigest(),
                coverage={"textual": "acquired-not-yet-attested-read",
                          "visual": "not-reviewed", "audio": "not-reviewed"})
    previous = old.get(vid, {})
    if previous.get("normalized_text_sha256") == meta[
            "normalized_text_sha256"]:
        meta["coverage"] = previous.get("coverage", meta["coverage"])
    return meta


def main():
    playlist = json.loads((PRIVATE / "playlist.json").read_text(
        encoding="utf-8-sig"))
    path = ROOT / "source-manifest.json"
    old = {}
    if path.exists():
        old = {r["video_id"]: r for r in json.loads(
            path.read_text(encoding="utf-8"))["recordings"]}
    entries = list(playlist["entries"]) + [SUPPLEMENT]
    results = [entry_record(i, e, old) for i, e in enumerate(entries, 1)]
    results[-1].update(role="official-supplement", playlist_position=None,
                       provenance=(
        "Official NPTEL course 111107164 outline lists this video as "
        "Lecture 44; it is absent from the supplied playlist."))
    manifest = {
        "course": "Advanced Linear Algebra", "institution": "IIT Roorkee",
        "instructor": "Prof. Premananda Bera",
        "playlist_url": PLAYLIST,
        "retrieved_at": datetime.now(timezone.utc).date().isoformat(),
        "expected_entries": len(results),
        "expected_teaching_recordings": len(results) - 1,
        "official_course_url": "https://nptel.ac.in/courses/111107164",
        "sequence_note": (
            "Playlist position 1 is the course promo. Positions 2-44 carry "
            "lectures 1-43. Lecture 44 is absent from the supplied playlist; "
            "the official NPTEL outline (course 111107164) lists it as video "
            "NYo2wbePAnI, appended here as position 61. Position 45 is Lecture 45; positions 46 and "
            "47 are Lectures 47 and 46 (uploaded out of order); positions "
            "48-53 are Lectures 48-53; positions 54 and 55 are Lectures 55 "
            "and 54 (out of order); positions 56-60 are Lectures 56-60."),
        "coverage_note": (
            "Caption acquisition is distinct from sequential reading. Chapter "
            "authors record reading scope in validation/source-read-*.json. "
            "No comprehensive audiovisual audit of boards or slides is "
            "claimed."),
        "recordings": results}
    path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False)
                    + "\n", encoding="utf-8", newline="\n")
    missing = [r["position"] for r in results if r["status"] != "acquired"]
    print(f"{len(results)} entries, missing captions: {missing}")


if __name__ == "__main__":
    main()
