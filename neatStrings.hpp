#ifndef NEAT_STRINGS
#define NEAT_STRINGS

#include <string>

/**
 * Contain custom toString()-esque and toWideString()-esque functions that mainly feature improved leading/trailing zero formatting
 */
namespace nstrings {
    /**
     * Get an arithmetic data type as a string
     * 
     * @returns The inputted data type formatted as a string
     */
    template <typename Type> std::string toString(Type input, bool specifyPositive = false) {
        // std::to_string() only works for arithmetic data types
        static_assert(std::is_arithmetic<Type>::value, "Type must be an arithmetic type");

        std::string output = std::to_string(input);
        // Remove any extra decimal places on the end
        if (std::is_floating_point<decltype(input)>::value) {
            for (unsigned long i = output.length() - 1; i >= 0; i--) {
                if (output[i] == '.') {
                    output.erase(i, 1);
                    break;
                } else if (output[i] != '0') {
                    break;
                }
                output.erase(i, 1);
            }
        }
        return (specifyPositive && input > 0 ? "+" : "") + output;
    }
    /**
     * Get an arithmetic data type as a string
     * This particular overload allows for extra formatting with leading/trailing zeros
     * 
     * @param beforeDecimal The amount of digits to place before the decimal point
     * @param afterDecimal The amount of digits to place after the decimal point (even for integral types)
     * @param add If true, beforeDecimal/afterDecimal switch from amount of digits to amount of leading/trailing zeros (regardless of digits present already)
     * @returns The inputted data type formatted as a string
     */
    template <typename Type> std::string toString_Places(Type input, unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false) {
        std::string output = str::toString(input, specifyPositive);
        unsigned long originalLength = output.length();

        if (beforeDecimal == 0 && afterDecimal == 0) {return output;}

        unsigned long beforeDecimalCount = 0;
        while (beforeDecimalCount < output.length()) {
            if (output[beforeDecimalCount] == '.') {break;}
            beforeDecimalCount++;
        }

        // Add leading zeros
        if (beforeDecimal > 0) {
            for (unsigned long i = add ? 0 : beforeDecimalCount; i < beforeDecimal; i++) {
                output.insert(input < 0 || specifyPositive ? 1 : 0, 1, '0');
            }

            if (afterDecimal == 0) {return output;}
        }

        unsigned long afterDecimalCount = output.find('.') == std::string::npos ? 0 : originalLength - beforeDecimalCount - 1;

        // Add trailing zeros
        if (afterDecimalCount == 0) {output += ".";}
        for (unsigned long i = add ? 0 : afterDecimalCount; i < afterDecimal; i++) {
            output += "0";
        }

        return output;
    }
    /**
     * Get an arithmetic data type as a string
     * This particular overload allows for extra formatting by specifying a length for each component
     * 
     * @param length The desired width of the output string
     * @param leading Whether to add the extra zeros to the beginning or end of the string
     * @returns The inputted data type formatted as a string
     */
    template <typename Type> std::string toString_Length(Type input, unsigned long length, bool leading = true, bool specifyPositive = false) {
        std::string output = str::toString(input, specifyPositive);
        bool hasDecimal = output.find('.') != std::string::npos;

        unsigned long beforeDecimal = 0;
        while (beforeDecimal < output.length()) {
            if (output[beforeDecimal] == '.') {break;}
            beforeDecimal++;
        }
        unsigned long afterDecimal = hasDecimal ? output.length() - beforeDecimal - 1 : 0;

        if (leading) {return str::toString_Places(input, length - afterDecimal - (hasDecimal ? 1 : 0), 0, false, specifyPositive);}
        return str::toString_Places(input, 0, length - beforeDecimal - 1, false, specifyPositive);
    }
    /**
     * Get a string as a wide string
     * 
     * @returns The inputted string as a wide string
     */
    std::wstring toWideString(std::string input) {
        std::wstring output;
        for (unsigned long i = 0; i < input.length(); i++) {
            output += input[i];
        }
        return output;
    }
    /**
     * Get an arithmetic data type as a wide string
     * 
     * @returns The inputted data type formatted as a wide string
     */
    template <typename Type> std::wstring toWideString(Type input, bool specifyPositive = false) {return str::toWideString(str::toString(input, specifyPositive));}
    /**
     * Get an arithmetic data type as a wide string
     * This particular overload allows for extra formatting with leading/trailing zeros
     * 
     * @param beforeDecimal The amount of digits to place before the decimal point
     * @param afterDecimal The amount of digits to place after the decimal point (even for integral types)
     * @param add If true, beforeDecimal/afterDecimal switch from amount of digits to amount of leading/trailing zeros (regardless of digits present already)
     * @returns The inputted data type formatted as a wide string
     */
    template <typename Type> std::wstring toWideString_Places(Type input, unsigned long beforeDecimal, unsigned long afterDecimal = 0, bool add = false, bool specifyPositive = false) {return str::toWideString(str::toString_Places(input, beforeDecimal, afterDecimal, add, specifyPositive));}
    /**
     * Get an arithmetic data type as a wide string
     * This particular overload allows for extra formatting by specifying a length for each component
     * 
     * @param length The desired width of the output wide string
     * @param leading Whether to add the extra zeros to the beginning or end of the wide string
     * @returns The inputted data type formatted as a wide string
     */
    template <typename Type> std::wstring toWideString_Length(Type input, unsigned long length, bool leading = true, bool specifyPositive = false) {return str::toWideString(str::toString_Length(input, length, leading, specifyPositive));}
}

#endif /* NEAT_STRINGS */
