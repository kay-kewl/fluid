#pragma once

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define COMMA ,
#define PAREN_OPEN (
#define PAREN_CLOSE )

#define CONCAT_IMPL(x, y) x ## y
#define CONCAT(x, y) CONCAT_IMPL(x, y)

