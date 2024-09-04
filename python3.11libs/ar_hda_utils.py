"""This module contains utility functions for easy HDA building and management,
including menu creation etc

NOTE: this module is intended for usage inside HDA sections, and as such requires 
the hou module available
"""

from typing import Union
import hou
import re
import itertools

RE_TOKEN_CAPITALIZE = r"\b(?P<token_name>\w{1,2})_?(?P<token_digit>\d)*\b"


def get_input_geo(node: hou.Node, input_idx: int = 0) -> hou.Geometry:
    """Given a node (usually from kwargs["node"]) returns the input geometry"""
    try:
        return node.input(input_idx).geometry()
    except hou.AttributeError as e:
        return None


def token_to_label(token: str) -> str:
    """Given a snake_case formatted token, converts it into a Title Case token."""

    token_slice = token.split("_")
    for i in range(len(token_slice)):
        token_slice[i] = (
            token_slice[i].upper()
            if len(token_slice[i]) <= 2
            else token_slice[i].capitalize()
        )

    return " ".join(token_slice)

    for m in re.finditer(RE_TOKEN_CAPITALIZE, token):
        token = (
            token[: m.start(0)]
            + f"{m['token_name'].upper()} {m['token_digit']}"
            + token[m.end(0) :]
        )

    return token.replace("_", " ")

    pieces = token.split("_")
    for p in pieces:
        m = re.search(RE_TOKEN_CAPITALIZE, p)
        if m is not None:
            print(m.groups())
            p = f"{m['token_name'].upper()} {m['token_digit']}"

    return " ".join(p)


def build_attrib_menu_listall(
    geo: hou.Geometry, attrib_type: hou.attribType
) -> list[str]:
    """Given an attribute type, returns a token, label entry list"""

    if geo is None:
        return []

    def _build_menu_list(attr_source: callable) -> list[str]:
        return [attr.name() for attr in attr_source() for _ in range(2)]

    if attrib_type == hou.attribType.Point:
        return _build_menu_list(geo.pointAttribs)

    if attrib_type == hou.attribType.Prim:
        return _build_menu_list(geo.primAttribs)

    if attrib_type == hou.attribType.Vertex:
        return _build_menu_list(geo.vertexAttribs)

    if attrib_type == hou.attribType.Global:
        return _build_menu_list(geo.globalAttribs)


def build_attrib_menu(
    geo: hou.Geometry,
    attrib_type: hou.attribType,
    attrib_data: Union[hou.attribData, list[hou.attribData]],
    size: Union[int, list[int]],
    include_array: bool = False,
) -> list[str]:
    """Given an attribute type, data and size, returns a [token, value] list"""

    if geo is None:
        return []

    def _build_menu_list(attr_source: callable) -> list[str]:
        return [
            attr.name()
            for attr in attr_source()
            for _ in range(2)
            if (
                (isinstance(attrib_data, list) and attr.dataType() in attrib_data)
                or attr.dataType() == attrib_data
            )
            and (
                (isinstance(size, list) and attr.size() in size) or attr.size() == size
            )
            and (attr.isArrayType() == include_array)
        ]

    if attrib_type == hou.attribType.Point:
        return _build_menu_list(geo.pointAttribs)

    if attrib_type == hou.attribType.Prim:
        return _build_menu_list(geo.primAttribs)

    if attrib_type == hou.attribType.Vertex:
        return _build_menu_list(geo.vertexAttribs)

    if attrib_type == hou.attribType.Global:
        return _build_menu_list(geo.globalAttribs)


def build_names_menu(geo: hou.Geometry) -> list[str]:
    """Returns a token, label entry list with all the unique named primitives"""

    if geo is None:
        return []

    return [name for name in geo.primStringAttribValues("name") for _ in range(2)]


def build_group_menu(
    geo: hou.Geometry,
    type: str,
) -> list[str]:
    """Given a string descriptor for a geometry type ('prim', 'point', 'edge' or 'vertex')
    returns a [name, name] list for menu dropdown
    """

    if geo is None:
        return []

    def _build_group_list(group_source: callable) -> list[str]:
        return [grp.name() for grp in group_source() for _ in range(2)]

    if type == "prim":
        return _build_group_list(geo.primGroups)

    if type == "point":
        return _build_group_list(geo.pointGroups)

    if type == "edge":
        return _build_group_list(geo.edgeGroups)

    if type == "vertex":
        return _build_group_list(geo.vertexGroups)


def build_entry_menu(entries: list[str], format_label: bool = False) -> list[str]:
    """Given a list of entry token for a menu, builds a [token, label] list"""

    if format_label:
        return list(
            itertools.chain.from_iterable(
                zip(entries, [token_to_label(e) for e in entries])
            )
        )
    else:
        return list(
            itertools.chain.from_iterable(itertools.repeat(e, 2) for e in entries)
        )


if __name__ == "__main__":
    tokens = ["uv", "uv2", "uv_2", "occlusion_amount"]
    print([token_to_label(t) for t in tokens])
