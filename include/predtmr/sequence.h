#pragma once

#include <string>
#include <vector>

namespace predtmr {

enum class SequenceFormat {
    unknown,
    plain,
    fasta,
    swissprot,
    pdb,
};

struct TransmembraneSegment {
    int start{0};
    int end{0};
    int extent{0};
    int direction{0};
};

class Sequence {
public:
    Sequence() = default;
    Sequence(std::string name, std::string data, SequenceFormat format)
        : name_(std::move(name)), data_(std::move(data)), format_(format) {}

    const std::string& name() const noexcept { return name_; }
    const std::string& data() const noexcept { return data_; }
    SequenceFormat format() const noexcept { return format_; }
    const std::vector<TransmembraneSegment>& observed_transmembrane() const noexcept {
        return observed_transmembrane_;
    }

    void set_name(std::string name) { name_ = std::move(name); }
    void set_data(std::string data) { data_ = std::move(data); }
    void set_format(SequenceFormat format) noexcept { format_ = format; }
    void add_observed_transmembrane(const TransmembraneSegment& segment) {
        observed_transmembrane_.push_back(segment);
    }

    bool is_transmembrane() const noexcept { return !observed_transmembrane_.empty(); }
    int molecule_length() const noexcept { return static_cast<int>(data_.size()); }

private:
    std::string name_;
    std::string data_;
    SequenceFormat format_{SequenceFormat::unknown};
    std::vector<TransmembraneSegment> observed_transmembrane_;
};

const char* format_name(SequenceFormat format) noexcept;

}  // namespace predtmr
