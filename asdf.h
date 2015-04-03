/*
 * Copyright (c) 2015 Claus Christmann <hcc |ä| gatech.edu>.  
 *   
 * Licensed under the Apache Lice*nse, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
