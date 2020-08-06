#include "TextTable.hpp"

std::ostream& operator<<(std::ostream& stream, TextTable const& table) {
  table.setup();
  stream << table.ruler() << "\n";
  for (auto rowIterator = table.rows().begin();
       rowIterator != table.rows().end();
       ++rowIterator) {
    TextTable::Row const& row = *rowIterator;
    stream << table.vertical();
    for (unsigned i = 0; i < row.size(); ++i) {
      auto alignment = table.alignment(i) == TextTable::Alignment::LEFT ? std::left : std::right;
      stream << std::setw(table.width(i)) << alignment << row[i];
      stream << table.vertical();
    }
    stream << "\n";
    stream << table.ruler() << "\n";
  }

  return stream;
}
