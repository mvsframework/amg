/*
 * Copyright (c) 2014, Claus Christmann <hcc@gatech.edu>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     1) Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *     2) Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 *
 */

/* inclusion guard */
#ifndef __ASDF_H__
#define __ASDF_H__

#include <glib-object.h>
/*
 * Potentially, include other headers on which this header depends.
 */

/*
 * Type macros.
 */
#define TYPE_ASDF                  (asdf_get_type ())
#define ASDF(obj)                  (G_TYPE_CHECK_INSTANCE_CAST ((obj), TYPE_ASDF, asdf))
#define IS_ASDF(obj)               (G_TYPE_CHECK_INSTANCE_TYPE ((obj), TYPE_ASDF))
#define ASDF_CLASS(klass)          (G_TYPE_CHECK_CLASS_CAST ((klass), TYPE_ASDF, asdfClass))
#define IS_ASDF_CLASS(klass)       (G_TYPE_CHECK_CLASS_TYPE ((klass), TYPE_ASDF))
#define ASDF_GET_CLASS(obj)        (G_TYPE_INSTANCE_GET_CLASS ((obj), TYPE_ASDF, asdfClass))

typedef struct _asdf       asdf;
typedef struct _asdfClass   asdfClass;

struct _asdf {
    GObject parent_instance;

    /* instance members */
};

struct _asdfClass {
    GObjectClass parent_class;
    /* class members */
};

/* used by TYPE_ASDF */
GType asdf_get_type ( void );

/*
 * Method definitions.
 */

#endif /* __ASDF_H__ */
