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
    // Every slot is a binary edge, so location[i] belongs to fragment[i]
    size_t n = fragment.size();

    std::vector<std::pair<Container, Container>>
        bodies_and_locations;

    bodies_and_locations.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        Container body;
        Container body_location;

        body.reserve(n - 1);

        body_location.reserve(n - 1);

        for (size_t j = i + 1;
             j < i + n;
             ++j)
        {
            size_t idx =
                j % n;

            body.push_back(
                fragment[idx]
            );

            body_location.push_back(
                location[idx]
            );
        }

        bodies_and_locations.emplace_back(
            std::move(body),
            std::move(body_location)
        );
    }

    return bodies_and_locations;
}



template<typename EdgeContainer1,
         typename EdgeContainer2>
std::vector<EdgeId>
direct_random_non_cyclic_selection(
    const std::vector<EdgeId>& edges,
    const EdgeContainer1& previous_edges,
    const EdgeContainer2& cyclic_binaries,
    size_t N
)
{
    std::vector<EdgeId> result;
    result.reserve(N);

    size_t attempts = 0;

    while (result.size() < N && attempts < 20)
    {
        size_t remaining =
            N - result.size();

        size_t batch_size =
            std::min(
                edges.size(),
                std::max<size_t>(
                    remaining * 2,
                    4
                )
            );

        auto candidates =
            get_random_elements(
                edges,
                batch_size
            );

        size_t before =
            result.size();

        for (EdgeId e : candidates)
        {
            if (has(previous_edges,e))
            {
                continue;
            }

            if (has(cyclic_binaries,e))
            {
                continue;
            }

            if (has(result,e))
            {
                continue;
            }

            result.push_back(e);

            if (result.size() == N)
            {
                break;
            }
        }

        ++attempts;
    }

    return result;
}




template<typename Container>
std::vector<Container>
find_canonical_orders(
    const Container& fragment)
{
    size_t n = fragment.size();

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
        for (size_t shift = 0; shift < n; ++shift) {
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
    // Every slot is a binary edge, so location[i] belongs to fragment[i]
    const size_t n = fragment.size();

    // Among equally canonical orders, pick the one with the lexicographically smallest location
    auto location_order_less =
        [&](const Container& lhs,
            const Container& rhs)
    {
        for (size_t i = 0; i < n; ++i) {
            if (location[lhs[i]] != location[rhs[i]]) {
                return location[lhs[i]] < location[rhs[i]];
            }
        }

        return false;
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

    Container ordered_location;

    ordered_location.resize(n);

    for (size_t i = 0; i < n; ++i) {
        ordered_fragment[i] = fragment[best_order[i]];
        ordered_location[i] = location[best_order[i]];
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
