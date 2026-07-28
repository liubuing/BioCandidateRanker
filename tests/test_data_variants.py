import pytest

from biocandidate.data.variants import apply_substitutions


def test_apply_substitutions_supports_wildtype_and_multiple_changes() -> None:
    reference = "ACDEFGHIK"

    assert apply_substitutions(reference, "WT") == reference
    assert apply_substitutions(reference, "A1V-F5Y") == "VCDEYGHIK"


@pytest.mark.parametrize(
    "notation, message",
    [
        ("A0V", "outside"),
        ("A10V", "outside"),
        ("C1V", "mismatch"),
        ("A1V-A1G", "duplicate"),
        ("A1", "invalid"),
        ("", "must not be empty"),
    ],
)
def test_apply_substitutions_rejects_ambiguous_variants(
    notation: str, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        apply_substitutions("ACDEFGHIK", notation)


def test_apply_substitutions_rejects_noncanonical_reference() -> None:
    with pytest.raises(ValueError, match="canonical"):
        apply_substitutions("ACDX", "WT")
