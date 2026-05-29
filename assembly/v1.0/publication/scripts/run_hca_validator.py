#!/usr/bin/env python3
"""HCA schema validator wrapper.

hca-schema-validator 0.5.0 has no CLI entry point (__main__.py missing).
This wrapper imports HCAValidator and runs it programmatically.

Usage:
    python3 run_hca_validator.py <file.h5ad> [--output <report.log>]
"""

import argparse
import sys
import traceback


def main():
    parser = argparse.ArgumentParser(description="Run HCA schema validator on an h5ad file")
    parser.add_argument("input", help="Path to h5ad file")
    parser.add_argument("--output", "-o", help="Path to write validation report (default: stdout)")
    args = parser.parse_args()

    try:
        from hca_schema_validator import HCAValidator
    except ImportError:
        print("ERROR: hca_schema_validator not installed", file=sys.stderr)
        sys.exit(1)

    print(f"HCA Validator: validating {args.input}")

    try:
        validator = HCAValidator()
        # validate_adata(h5ad_path) -> bool; errors on validator.errors list
        is_valid = validator.validate_adata(args.input)
    except Exception:
        msg = f"ERROR: HCA validation failed with exception:\n{traceback.format_exc()}"
        if args.output:
            with open(args.output, "w") as f:
                f.write(msg)
        print(msg)
        sys.exit(1)

    # Format output
    lines = [f"is_valid: {is_valid}"]

    errors = getattr(validator, "errors", [])
    warnings = getattr(validator, "warnings", [])

    if errors:
        lines.append(f"errors ({len(errors)}):")
        for err in errors:
            lines.append(f"  ERROR: {err}")
    if warnings:
        lines.append(f"warnings ({len(warnings)}):")
        for warn in warnings:
            lines.append(f"  WARNING: {warn}")

    if not errors and not warnings:
        lines.append("No errors or warnings.")

    output_text = "\n".join(lines)

    if args.output:
        with open(args.output, "w") as f:
            f.write(output_text + "\n")

    print(output_text)
    sys.exit(0 if is_valid else 1)


if __name__ == "__main__":
    main()
