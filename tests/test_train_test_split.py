"""Tests for train/test split logic."""

import random

from train_test_split import (
    build_tree_from_branches,
    get_all_tips_iterative,
    iterative_test_selection,
)


def make_branches_dict(parent_of, hamming_of):
    """Convert fixture dicts to the branches_dict format
    expected by train_test_split: child -> (parent, hamming, train_test)."""
    branches = {}
    for child, parent in parent_of.items():
        hamming = hamming_of.get((parent, child), 0)
        branches[child] = (parent, hamming, "")
    return branches


class TestBuildTree:
    """Test tree construction from branches."""

    def test_basic_tree(self, parent_of, hamming_of):
        branches = make_branches_dict(parent_of, hamming_of)
        children, parents, root = build_tree_from_branches(branches)
        assert root == "X"
        assert set(parents.keys()) == {"Y", "A", "B", "C"}
        assert "Y" in children["X"]
        assert "C" in children["X"]
        assert "A" in children["Y"]
        assert "B" in children["Y"]

    def test_finds_correct_root(self, parent_of, hamming_of):
        branches = make_branches_dict(parent_of, hamming_of)
        _, _, root = build_tree_from_branches(branches)
        assert root == "X"


class TestGetAllTips:
    """Test tip discovery."""

    def test_finds_all_tips(self, parent_of, hamming_of):
        branches = make_branches_dict(parent_of, hamming_of)
        children, _, root = build_tree_from_branches(branches)
        tips = get_all_tips_iterative(children, root)
        assert set(tips) == {"A", "B", "C"}

    def test_no_internal_nodes(self, parent_of, hamming_of):
        branches = make_branches_dict(parent_of, hamming_of)
        children, _, root = build_tree_from_branches(branches)
        tips = get_all_tips_iterative(children, root)
        assert "X" not in tips
        assert "Y" not in tips


class TestIterativeTestSelection:
    """Test the test selection algorithm on a synthetic tree."""

    def _make_balanced_tree(self, n_tips=100):
        """Create a balanced binary tree with approximately n_tips leaves."""
        branches = {}
        node_id = 0
        root = f"node_{node_id}"
        current_layer = [root]

        while True:
            next_layer = []
            for parent in current_layer:
                for _ in range(2):
                    node_id += 1
                    child = f"node_{node_id}"
                    branches[child] = (parent, 2, "")
                    next_layer.append(child)
            # Stop when the next layer would be leaves with enough tips
            if len(next_layer) >= n_tips:
                break
            current_layer = next_layer

        return branches

    def test_all_tips_assigned(self):
        """Every tip is either train or test, with no overlap."""
        branches = self._make_balanced_tree(64)
        children, parents, root = build_tree_from_branches(branches)
        rng = random.Random(42)
        test_nodes, test_tips, total_tips = iterative_test_selection(
            children, parents, root, branches,
            test_proportion=0.2,
            mutations_back=3,
            max_clade_proportion=0.1,
            rng=rng,
        )

        all_tips = set(get_all_tips_iterative(children, root))
        train_tips = all_tips - test_tips
        # Every tip is either train or test
        assert train_tips | test_tips == all_tips
        assert len(train_tips & test_tips) == 0

    def test_produces_test_tips(self):
        """Selection actually marks some tips as test."""
        branches = self._make_balanced_tree(64)
        children, parents, root = build_tree_from_branches(branches)
        rng = random.Random(42)
        _, test_tips, total_tips = iterative_test_selection(
            children, parents, root, branches,
            test_proportion=0.2,
            mutations_back=3,
            max_clade_proportion=0.1,
            rng=rng,
        )
        assert len(test_tips) > 0
        assert len(test_tips) <= total_tips

    def test_deterministic_with_seed(self):
        """Same seed produces the same split."""
        branches = self._make_balanced_tree(64)
        children, parents, root = build_tree_from_branches(branches)

        rng1 = random.Random(42)
        _, test_tips_1, _ = iterative_test_selection(
            children, parents, root, branches,
            test_proportion=0.2, mutations_back=3,
            max_clade_proportion=0.1, rng=rng1,
        )

        rng2 = random.Random(42)
        _, test_tips_2, _ = iterative_test_selection(
            children, parents, root, branches,
            test_proportion=0.2, mutations_back=3,
            max_clade_proportion=0.1, rng=rng2,
        )

        assert test_tips_1 == test_tips_2
