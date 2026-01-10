#!/usr/bin/env python3
"""Extract color mappings from Auspice JSON for use in visualizations."""

import argparse
import json


def extract_colors(auspice_json):
    """
    Extract color mappings from Auspice JSON metadata.

    Returns a dict with structure:
    {
        "field_key": {
            "title": "Display Title",
            "colors": {"value1": "#hex1", "value2": "#hex2", ...}
        }
    }
    """
    colors_data = {}

    # Check for meta.colorings
    if 'meta' not in auspice_json or 'colorings' not in auspice_json['meta']:
        return colors_data

    # Process each coloring
    for coloring in auspice_json['meta']['colorings']:
        # Only process categorical colorings with scales
        if coloring.get('type') == 'categorical' and 'scale' in coloring:
            key = coloring.get('key')
            title = coloring.get('title', key)

            # Convert scale array to dict
            color_map = {}
            for value, color in coloring['scale']:
                color_map[value] = color

            colors_data[key] = {
                'title': title,
                'colors': color_map
            }

    return colors_data


def main():
    parser = argparse.ArgumentParser(description="Extract color mappings from Auspice JSON")
    parser.add_argument("--json", required=True, help="Path to the auspice.json file")
    parser.add_argument("--output", required=True, help="Output colors.json file")

    args = parser.parse_args()

    # Load Auspice JSON
    with open(args.json, 'r') as f:
        auspice_data = json.load(f)

    # Extract colors
    colors_data = extract_colors(auspice_data)

    # Write output
    with open(args.output, 'w') as f:
        json.dump(colors_data, f, indent=2)

    # Print summary
    if colors_data:
        print(f"Extracted color mappings for {len(colors_data)} fields:")
        for key, data in colors_data.items():
            print(f"  - {key}: {data['title']} ({len(data['colors'])} values)")
    else:
        print("No categorical color mappings found in Auspice JSON")


if __name__ == "__main__":
    main()