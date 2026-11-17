"""Build a QGIS plugin zip and a plugins.xml for a custom plugin repository.

Used by .github/workflows/qgis-plugin-release.yml when a GitHub release is
published. Files marked export-ignore in .gitattributes are left out of the zip.

    python .github/scripts/build_qgis_release.py --plugin-dir NAME --tag v1.2.3 \
        --repo owner/repo --out dist
"""
import argparse
import configparser
import datetime
import os
import subprocess
import sys
from xml.sax.saxutils import escape, quoteattr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plugin-dir", required=True,
                    help="folder name the plugin is installed under (the plugin key)")
    ap.add_argument("--tag", required=True, help="release tag, e.g. v0.2.24")
    ap.add_argument("--repo", required=True, help="owner/repo on GitHub")
    ap.add_argument("--out", default="dist")
    args = ap.parse_args()

    cfg = configparser.ConfigParser(interpolation=None)
    cfg.read("metadata.txt", encoding="utf-8")
    md = cfg["general"]
    version = md["version"].strip()

    if args.tag.lstrip("vV") != version:
        sys.exit(
            "Release tag %s does not match version=%s in metadata.txt. "
            "Update metadata.txt or use tag v%s." % (args.tag, version, version)
        )

    os.makedirs(args.out, exist_ok=True)
    zip_name = args.plugin_dir + ".zip"
    zip_path = os.path.join(args.out, zip_name)
    subprocess.run(
        ["git", "archive", "--format=zip", "--prefix=%s/" % args.plugin_dir,
         "-o", zip_path, "HEAD"],
        check=True,
    )

    def get(key, default=""):
        return md.get(key, default).strip()

    url = "https://github.com/%s/releases/download/%s/%s" % (args.repo, args.tag, zip_name)
    today = datetime.date.today().isoformat()
    fields = [
        ("description", get("description")),
        ("about", get("about")),
        ("version", version),
        ("qgis_minimum_version", get("qgisMinimumVersion", "3.0")),
        ("qgis_maximum_version", get("qgisMaximumVersion")),
        ("homepage", get("homepage")),
        ("file_name", zip_name),
        ("icon", get("icon")),
        ("author_name", get("author")),
        ("download_url", url),
        ("uploaded_by", get("author")),
        ("create_date", today),
        ("update_date", today),
        ("experimental", get("experimental", "False")),
        ("deprecated", get("deprecated", "False")),
        ("tracker", get("tracker")),
        ("repository", get("repository")),
        ("tags", get("tags")),
        ("server", get("server", "False")),
    ]
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        "<plugins>",
        "  <pyqgis_plugin name=%s version=%s plugin_id=%s>"
        % (quoteattr(get("name", args.plugin_dir)), quoteattr(version), quoteattr(args.plugin_dir)),
    ]
    for tag, value in fields:
        if value:
            lines.append("    <%s>%s</%s>" % (tag, escape(value), tag))
    lines += ["  </pyqgis_plugin>", "</plugins>", ""]
    with open(os.path.join(args.out, "plugins.xml"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines))
    print("Built %s and plugins.xml for %s %s" % (zip_path, args.plugin_dir, version))


if __name__ == "__main__":
    main()
