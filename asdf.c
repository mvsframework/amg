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

#include "asdf.h"

/*
 * forward definitions
 */
G_DEFINE_TYPE ( asdf, asdf, G_TYPE_OBJECT );
/*
/* forward declarations of default virtual methods
 */

static void
asdf_dispose ( GObject *gobject )
{
    asdf *self = ASDF ( gobject );

    /*
     * In dispose, you are supposed to free all types referenced from this
     * object which might themselves hold a reference to self. Generally,
     * the most simple solution is to unref all members on which you own a
     * reference.
     */

    /* Chain up to the parent class */
    G_OBJECT_CLASS ( asdf_parent_class )->dispose ( gobject );
}

static void
asdf_finalize ( GObject *gobject )
{
    asdf *self = ASDF ( gobject );

    /* Chain up to the parent class */
    G_OBJECT_CLASS ( asdf_parent_class )->finalize ( gobject );
}

static void
asdf_init ( asdf *self )
{
    /* initialize all public and private members to reasonable default values. */

    /* Default implementations for virtual methods
     * For pure-virtual functions, set these to NULL
     */
}
