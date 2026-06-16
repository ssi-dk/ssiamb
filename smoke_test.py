from __future__ import annotations

from ssiamb.cli import main


def test_smoke() -> None:
    try:
        main(["--version"])
    except SystemExit as exc:
        assert exc.code == 0


if __name__ == "__main__":
    test_smoke()
