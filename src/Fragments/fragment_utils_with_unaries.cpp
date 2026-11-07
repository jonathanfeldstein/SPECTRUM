#ifndef SPECTRUM_FRAGMENT_UTILS_H
#define SPECTRUM_FRAGMENT_UTILS_H

// #include "../Utils/qrand.h"

template<class T>
std::vector<T>
get_random_elements(
    const std::vector<T>& input,
    size_t N)
{
    thread_local std::mt19937 gen(
        std::random_device{}()
    );

    const size_t M = input.size();

    if (N >= M)
    {
        return input;
    }

    auto rand_bounded =
        [&](size_t bound) -> size_t
    {
        return gen() % bound;
    };

    // -------------------------------------------------
    // N == 1
    // -------------------------------------------------

    if (N == 1)
    {
        return {
            input[rand_bounded(M)]
        };
    }

    // -------------------------------------------------
    // N == 2
    // -------------------------------------------------

    if (N == 2)
    {
        size_t a = rand_bounded(M);

        size_t b = rand_bounded(M - 1);

        if (b >= a)
        {
            ++b;
        }

        return {
            input[a],
            input[b]
        };
    }

    // -------------------------------------------------
    // N == 3
    // -------------------------------------------------

    if (N == 3)
    {
        size_t a = rand_bounded(M);

        size_t b = rand_bounded(M - 1);

        if (b >= a)
        {
            ++b;
        }

        size_t c = rand_bounded(M - 2);

        size_t lo = std::min(a, b);
        size_t hi = std::max(a, b);

        if (c >= lo)
        {
            ++c;
        }

        if (c >= hi)
        {
            ++c;
        }

        return {
            input[a],
            input[b],
            input[c]
        };
    }

    // -------------------------------------------------
    // N == 4
    // -------------------------------------------------

    if (N == 4)
    {
        std::array<size_t,4> idx;

        idx[0] = rand_bounded(M);

        idx[1] = rand_bounded(M - 1);
        if (idx[1] >= idx[0])
        {
            ++idx[1];
        }

        idx[2] = rand_bounded(M - 2);

        {
            std::array<size_t,2> used = {
                idx[0],
                idx[1]
            };

            std::sort(used.begin(), used.end());

            if (idx[2] >= used[0]) ++idx[2];
            if (idx[2] >= used[1]) ++idx[2];
        }

        idx[3] = rand_bounded(M - 3);

        {
            std::array<size_t,3> used = {
                idx[0],
                idx[1],
                idx[2]
            };

            std::sort(used.begin(), used.end());

            if (idx[3] >= used[0]) ++idx[3];
            if (idx[3] >= used[1]) ++idx[3];
            if (idx[3] >= used[2]) ++idx[3];
        }

        return {
            input[idx[0]],
            input[idx[1]],
            input[idx[2]],
            input[idx[3]]
        };
    }

    // -------------------------------------------------
    // M - N == 1
    // -------------------------------------------------

    if (M - N == 1)
    {
        size_t skip = rand_bounded(M);

        std::vector<T> result;
        result.reserve(N);

        for (size_t i = 0; i < M; ++i)
        {
            if (i != skip)
            {
                result.push_back(input[i]);
            }
        }

        return result;
    }

    // -------------------------------------------------
    // M - N == 2
    // -------------------------------------------------

    if (M - N == 2)
    {
        size_t a = rand_bounded(M);

        size_t b = rand_bounded(M - 1);

        if (b >= a)
        {
            ++b;
        }

        std::vector<T> result;
        result.reserve(N);

        for (size_t i = 0; i < M; ++i)
        {
            if (i != a && i != b)
            {
                result.push_back(input[i]);
            }
        }

        return result;
    }

    thread_local std::vector<size_t> indices;

    indices.resize(M);

    std::iota(
        indices.begin(),
        indices.end(),
        0
    );

    for (size_t i = 0; i < N; ++i)
    {
        size_t j =
            i +
            rand_bounded(M - i);

        std::swap(
            indices[i],
            indices[j]
        );
    }

    std::vector<T> result;
    result.reserve(N);

    for (size_t i = 0; i < N; ++i)
    {
        result.push_back(
            input[indices[i]]
        );
    }

    return result;
 

}


template<typename Container>
Container order_fragment(const Container& fragment)
{
    Container reversed_fragment = fragment;

    std::reverse(
        reversed_fragment.begin(),
        reversed_fragment.end()
    );

    if (fragment < reversed_fragment) {
        return fragment;
    }

    return reversed_fragment;
}
template<typename Container>
std::vector<std::pair<Container, Container>>
find_bodies_of_cyclic_fragment(
    const Container& fragment,
    const Container& location)
{
    size_t n = fragment.size();

    if (n % 2 != 0) {
        return {};
    }

    std::vector<std::pair<Container, Container>>
        bodies_and_locations;

    bodies_and_locations.reserve(n);

    std::vector<int> location_template(
        n,
        -1
    );

    size_t k = 0;

    for (size_t j = 0; j < n; ++j)
    {
        if (fragment[j] != 0)
        {
            location_template[j] =
                static_cast<int>(location[k]);

            ++k;
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        Container body;
        Container body_location;

        body.reserve(n + 1);

        body_location.reserve(
            location.size()
        );

        for (size_t j = i + 1;
             j < i + n;
             ++j)
        {
            size_t idx =
                j % n;

            body.push_back(
                fragment[idx]
            );

            if (location_template[idx] != -1)
            {
                body_location.push_back(
                    static_cast<typename Container::value_type>(
                        location_template[idx]
                    )
                );
            }
        }

        if (i % 2 == 0)
        {
            body.insert(
                body.begin(),
                0
            );

            body.push_back(0);
        }

        bodies_and_locations.emplace_back(
            std::move(body),
            std::move(body_location)
        );
    }

    return bodies_and_locations;
}

template<typename Container>
std::vector<Container>
find_canonical_orders(
    const Container& fragment)
{
    size_t n = fragment.size();

    if (n % 2 != 0) {
        return {};
    }

    std::vector<Container> result;

    auto compare_less =
        [&](const Container& a,
            const Container& b)
    {
        for (size_t i = 0; i < n; ++i) {
            if (fragment[a[i]] != fragment[b[i]])
            {
                return
                    fragment[a[i]]
                    <
                    fragment[b[i]];
            }
        }

        return false;
    };

    auto compare_equal =
        [&](const Container& a,
            const Container& b)
    {
        for (size_t i = 0; i < n; ++i) {
            if (fragment[a[i]] != fragment[b[i]]) {
                return false;
            }
        }

        return true;
    };

    Container canonical;

    bool first = true;

    for (int reflect = 0; reflect < 2; ++reflect) {
        for (size_t shift = 0; shift < n; shift += 2) {
            Container permutation;

            permutation.resize(n);

            for (size_t i = 0; i < n; ++i) {
                size_t idx = (i + shift) % n;

                if (reflect) {
                    idx = (n - idx) % n;
                }

                permutation[i] = idx;
            }

            if (first || compare_less(permutation, canonical)) {
                canonical = permutation;

                result = {permutation};

                first = false;
            } else if (compare_equal(permutation, canonical)) {
                result.push_back(permutation);
            }
        }
    }

    return result;
}

template<typename Container>
std::pair<Container, Container>
reorder_cyclic_fragment(
    const Container& fragment,
    const Container& location,
    const std::vector<Container>& possible_orders)
{
    const size_t n = fragment.size();

    if (n % 2 != 0) {
        return {};
    }

    std::vector<int> location_template(n,-1);

    size_t k = 0;

    for (size_t i = 0; i < n; ++i) {
        if (fragment[i] != 0) {
            location_template[i] = static_cast<int>(location[k]);
            ++k;
        }
    }

    auto location_order_less =
        [&](const Container& lhs,
            const Container& rhs)
    {
        size_t lhs_pos = 0;
        size_t rhs_pos = 0;

        while (true)
        {
            while (lhs_pos < n && location_template[lhs[lhs_pos]] == -1) {
                ++lhs_pos;
            }

            while (rhs_pos < n && location_template[rhs[rhs_pos]] == -1) {
                ++rhs_pos;
            }

            bool lhs_end = lhs_pos == n;

            bool rhs_end = rhs_pos == n;

            if (lhs_end && rhs_end) {
                return false;
            }

            int lhs_value = location_template[lhs[lhs_pos]];

            int rhs_value = location_template[rhs[rhs_pos]];

            if (lhs_value < rhs_value) {
                return true;
            }

            if (lhs_value > rhs_value) {
                return false;
            }

            ++lhs_pos;
            ++rhs_pos;
        }
    };

    size_t best_order_idx = 0;

    for (size_t i = 1; i < possible_orders.size(); ++i) {
        if (location_order_less(possible_orders[i], possible_orders[best_order_idx])) {
            best_order_idx = i;
        }
    }

    const auto& best_order = possible_orders[best_order_idx];

    Container ordered_fragment;

    ordered_fragment.resize(n);

    for (size_t i = 0; i < n; ++i) {
        ordered_fragment[i] = fragment[best_order[i]];
    }

    Container ordered_location;

    ordered_location.reserve(location.size());

    for (size_t i = 0; i < n; ++i) {
        int loc = location_template[best_order[i]];

        if (loc != -1) {
            ordered_location.push_back(static_cast<typename Container::value_type>(loc));
        }
    }

    return {std::move(ordered_fragment), std::move(ordered_location)};
}

template<typename Container>
std::pair<Container, Container>
order_cyclic_fragment_and_location(
    Container fragment,
    Container location)
{
    auto ordered_indices_possibilities = find_canonical_orders(fragment);

    return reorder_cyclic_fragment(fragment, location, ordered_indices_possibilities);
}

#endif //SPECTRUM_FRAGMENT_UTILS_H
