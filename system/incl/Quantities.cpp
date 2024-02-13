/*
 * Static function to find quantity by name
 * Author: JJ
 * Date 2022/09/09
 */
#include <algorithm>
#include "Quantities.hpp"

int quant::idx_from_code(std::string word)
{
    std::vector<std::string>::const_iterator quantit;
    quantit = std::find(quant::code.begin(), quant::code.end(), word);
    if (quantit == quant::code.end())
    {
        return -1;
    }
    return std::distance(quant::code.begin(), quantit);
}